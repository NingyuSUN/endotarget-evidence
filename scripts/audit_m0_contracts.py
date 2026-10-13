#!/usr/bin/env python3
"""Check real donor count/design interfaces; this alone cannot pass scientific M0."""
import csv
import datetime
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np

class ContractError(ValueError):
    pass

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def table(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def clinical_contract(rows, minimum_per_group, cell_field):
    if not rows:
        raise ContractError("empty clinical ledger")
    ids = [r["donor_id"] for r in rows]
    if len(ids) != len(set(ids)):
        raise ContractError("duplicate biological donor in clinical ledger")
    groups = {g: sum(r["disease"] == g for r in rows) for g in ("ENDO", "non-ENDO")}
    if sum(groups.values()) != len(rows) or min(groups.values()) < minimum_per_group:
        raise ContractError("unknown disease category or insufficient eligible donors")
    if any(int(r[cell_field]) < 50 for r in rows):
        raise ContractError("primary cell threshold failed")
    return ids, groups

def count_contract(path, donor_ids, integer_required):
    totals = [Decimal(0) for _ in donor_ids]
    gene_ids = []
    fractional_values = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, [])
        if header != ["gene_id"] + donor_ids:
            raise ContractError("count header/order differs from the frozen clinical donor order")
        seen = set()
        for row in reader:
            if len(row) != len(header):
                raise ContractError("count row width mismatch")
            gene = row[0]
            if not re.fullmatch(r"ENSG[0-9]{11}", gene) or gene in seen:
                raise ContractError("unstable or duplicate gene ID")
            seen.add(gene)
            gene_ids.append(gene)
            for j, token in enumerate(row[1:]):
                try:
                    value = Decimal(token)
                except InvalidOperation as e:
                    raise ContractError("unparseable count") from e
                if not value.is_finite() or value < 0:
                    raise ContractError("negative or nonfinite count")
                if value != value.to_integral_value():
                    fractional_values += 1
                    if integer_required:
                        raise ContractError("fractional value in declared integer HECA input")
                totals[j] += value
    if not gene_ids or any(x <= 0 for x in totals):
        raise ContractError("empty gene matrix or zero-total donor")
    return {
        "genes": len(gene_ids), "donors": len(donor_ids),
        "donor_order": donor_ids, "fractional_entries": fractional_values,
        "library_totals": {d: str(t) for d, t in zip(donor_ids, totals)},
        "integer_required": integer_required,
    }, set(gene_ids)

def design_contract(path, clinical_rows):
    with Path(path).open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)
    ids = [r["donor_id"] for r in clinical_rows]
    if not fields or fields[0] != "donor_id" or [r["donor_id"] for r in rows] != ids:
        raise ContractError("design donor order mismatch")
    columns = fields[1:]
    if "ENDO_minus_non_ENDO" not in columns or len(columns) != len(set(columns)):
        raise ContractError("missing disease coefficient or duplicate design column")
    try:
        x = np.asarray([[float(r[c]) for c in columns] for r in rows])
    except ValueError as e:
        raise ContractError("nonnumeric design value") from e
    if not np.isfinite(x).all():
        raise ContractError("nonfinite design")
    expected = np.asarray([int(r["disease"] == "ENDO") for r in clinical_rows])
    if not np.array_equal(x[:, columns.index("ENDO_minus_non_ENDO")], expected):
        raise ContractError("design disease labels disagree with clinical ledger")
    rank = int(np.linalg.matrix_rank(x))
    df = len(rows) - rank
    if rank != x.shape[1] or df < 5:
        raise ContractError("rank deficiency or residual degrees of freedom below five")
    return {"rank": rank, "columns": columns, "residual_df": df, "donors": len(rows)}

def verify_receipt_output(receipt_path, counts_path):
    receipt = json.loads(Path(receipt_path).read_text())
    if receipt.get("normalized") is not False:
        raise ContractError("count receipt must explicitly declare normalized=false")
    records = [r for r in receipt.get("outputs", []) if r.get("path") == str(counts_path)]
    if len(records) != 1:
        raise ContractError("count output missing/duplicated in provenance receipt")
    record = records[0]
    if sha256(counts_path) != record["sha256"] or Path(counts_path).stat().st_size != record["bytes"]:
        raise ContractError("count file differs from preparation receipt")
    for path, expected in receipt.get("input_hashes", {}).items():
        if sha256(path) != expected:
            raise ContractError("upstream input changed since preparation: " + path)
    return receipt

def main():
    specs = {
        "discovery": (8, "primary_lineage_cells"),
        "heca": (5, "primary_nuclei"),
    }
    result = {
        "checked_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "status": "COUNT_DESIGN_CONTRACTS_INCOMPLETE",
        "scientific_m0_passed": False,
        "scientific_effects_computed": False,
        "cohorts": {}, "failures": [],
        "scope": "Structural count and clinical-design interfaces, not source biological validation or full M0",
    }
    genes = {}
    frozen_paths = [Path(__file__)]
    for cohort, (minimum, cell_field) in specs.items():
        try:
            clinical = Path("manifests/final_" + cohort + "_clinical_selected.tsv")
            counts = Path("data/processed") / cohort / "counts_stable.tsv.gz"
            design = Path("reports/final_clinical_design_" + cohort + "_primary.tsv")
            receipt_path = Path("receipts/" + cohort + "_counts_preparation.json")
            rows = table(clinical)
            ids, groups = clinical_contract(rows, minimum, cell_field)
            receipt = verify_receipt_output(receipt_path, counts)
            values, genes[cohort] = count_contract(counts, ids, cohort == "heca")
            design_result = design_contract(design, rows)
            result["cohorts"][cohort] = {
                "counts": values, "groups": groups, "design": design_result,
                "status": "COUNT_DESIGN_CONTRACT_PASS",
                "normalization": "none",
                "clinical_path": str(clinical), "counts_path": str(counts),
            }
            frozen_paths += [clinical, counts, design, receipt_path]
            frozen_paths += [Path(x) for x in receipt.get("input_hashes", {})]
        except (OSError, ValueError, KeyError) as e:
            result["failures"].append({"cohort": cohort, "error": type(e).__name__ + ": " + str(e)})
    if not result["failures"]:
        result["status"] = "COUNT_DESIGN_CONTRACTS_PASS_NOT_FULL_M0"
        result["label_blind_gene_coverage"] = {
            "shared_stable_ids": len(genes["discovery"] & genes["heca"]),
            "discovery_only_ids": len(genes["discovery"] - genes["heca"]),
            "heca_only_ids": len(genes["heca"] - genes["discovery"]),
            "interpretation": "ID availability only; not expression filtering or the comparative graph universe",
        }
    root = Path(".").resolve()
    locked = []
    for path in sorted({p.resolve() for p in frozen_paths}):
        try:
            relative = path.relative_to(root)
        except ValueError as e:
            raise ContractError("input path escaped the authorized project root") from e
        if path.is_file():
            locked.append({"path": str(relative), "bytes": path.stat().st_size, "sha256": sha256(path)})
    result["verified_inputs"] = locked
    result["runtime"] = {"python": sys.version.split()[0], "numpy": np.__version__}
    Path("receipts/m0_count_design_contracts.json").write_text(json.dumps(result, indent=2) + "\n")
    with Path("manifests/m0_count_design_inputs.lock.tsv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["path", "bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        writer.writerows(locked)
    print(json.dumps({k:v for k,v in result.items() if k != "verified_inputs" and k != "cohorts"}, indent=2))
    return 1 if result["failures"] else 0

if __name__ == "__main__":
    raise SystemExit(main())
