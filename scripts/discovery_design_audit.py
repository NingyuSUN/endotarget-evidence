#!/usr/bin/env python3
"""Metadata-only design checks; never reads expression or fits disease effects."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
import numpy as np

INPUT = Path("manifests/discovery_clinical_candidates.tsv")
COVERAGE = Path("data/input_routes/osf/independent_donor_coverage.csv")


def design(rows, terms):
    columns = ["intercept", "ENDO_minus_non_ENDO"]
    values = [np.ones(len(rows)), np.array([r["endometriosis status"] == "ENDO" for r in rows], dtype=float)]
    for term in terms:
        if term == "age":
            a = np.array([float(r["age_years"]) for r in rows])
            values.append(a - a.mean())
            columns.append("age_centered_years")
        elif term == "adenomyosis":
            values.append(np.array([int(r["Adenomyosis"]) for r in rows], dtype=float))
            columns.append("adenomyosis")
        else:
            field = {"capture": "10x capturing batch", "sequencing": "sequencing batch"}[term]
            levels = sorted({r[field] for r in rows})
            if not all(levels):
                raise ValueError("missing batch metadata")
            for level in levels[1:]:
                columns.append(term + ":" + level + "_vs_" + levels[0])
                values.append(np.array([r[field] == level for r in rows], dtype=float))
    x = np.column_stack(values)
    if not np.isfinite(x).all():
        raise ValueError("nonfinite design")
    return x, columns


def summarize(rows, terms):
    x, columns = design(rows, terms)
    rank = int(np.linalg.matrix_rank(x))
    nuisance = np.delete(x, 1, axis=1)
    disease_residual = x[:, 1] - nuisance @ np.linalg.lstsq(nuisance, x[:, 1], rcond=None)[0]
    counts = Counter(r["endometriosis status"] for r in rows)
    contrast_identifiable = np.linalg.matrix_rank(x) > np.linalg.matrix_rank(nuisance)
    return {"terms": list(terms), "n": len(rows), "groups": dict(counts), "columns": columns,
            "rank": rank, "n_columns": int(x.shape[1]), "full_rank": rank == x.shape[1],
            "residual_df": len(rows) - rank, "disease_contrast_identifiable": bool(contrast_identifiable),
            "disease_information_unweighted": float(disease_residual @ disease_residual),
            "count_minimum_pass": min(counts.get("ENDO", 0), counts.get("non-ENDO", 0)) >= 8,
            "structural_gate_pass": bool(rank == x.shape[1] and len(rows) - rank >= 5 and contrast_identifiable and min(counts.get("ENDO", 0), counts.get("non-ENDO", 0)) >= 8),
            "scope": "metadata-only algebra; not scientific validation, not frozen model"}


def main():
    rows = list(csv.DictReader(INPUT.open(), delimiter="\t"))
    coverage = {r["sample"]: r for r in csv.DictReader(COVERAGE.open()) if r["population_id"] == "P077"}
    strict = [r for r in rows if r["proposal_conservative_clinical_candidate"] == "True"]
    primary = [r for r in rows if r["author_proliferative_candidate"] == "True" and r["Minor exclusion criteria logical"] == "0" and not r["clinical_disease_inconsistency"]]
    assert len(primary) == 22 and len({r["sample_id"] for r in primary}) == 22
    assert len(strict) == 20
    for r in primary:
        assert r["sample_id"] in coverage
        assert coverage[r["sample_id"]]["GSM"] == r["GSM"]
        assert coverage[r["sample_id"]]["disease"] == r["endometriosis status"]
        assert int(coverage[r["sample_id"]]["cells_independent"]) >= 50
    scenarios = {"author_adjudicated_phase22": primary,
                 "strict_histology_concordance20": strict,
                 "author_phase_adenomyosis_free18": [r for r in primary if r["Adenomyosis"] == "0"],
                 "strict_phase_adenomyosis_free16": [r for r in strict if r["Adenomyosis"] == "0"]}
    models = [["age", "capture"], ["age", "capture", "adenomyosis"],
              ["age", "sequencing"], ["age", "sequencing", "adenomyosis"],
              ["age", "capture", "sequencing"], ["age"]]
    results = {k: [summarize(v, m) for m in models] for k, v in scenarios.items()}
    # Meaningful guard: nuisance perfectly equal to disease makes that contrast unidentifiable.
    test = [dict(primary[i], **{"10x capturing batch": primary[i]["endometriosis status"]}) for i in range(len(primary))]
    assert not summarize(test, ["capture"])["disease_contrast_identifiable"]
    assert all(not item["count_minimum_pass"] for item in results["strict_phase_adenomyosis_free16"])
    batches = {}
    for field in ["10x capturing batch", "sequencing batch"]:
        batches[field] = {b: dict(Counter(r["endometriosis status"] for r in primary if r[field] == b)) for b in sorted({r[field] for r in primary})}
    report = {"status": "PROPOSED_INPUTS_NOT_FROZEN", "numpy_version": np.__version__,
              "source_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [INPUT, COVERAGE]},
              "selected_ids": [r["sample_id"] for r in primary], "models": results, "batch_overlap": batches,
              "coverage": {"cells_min": min(int(coverage[r["sample_id"]]["cells_independent"]) for r in primary), "all_22_meet_50": True},
              "validation": {"GSM_disease_coverage_join": "PASS", "perfect_confounding_guard": "PASS", "minimum_donor_guard": "PASS"},
              "interpretation": "Rank and residual degrees of freedom assess estimability, not biological validity, power, or removal of residual confounding. Sequencing batch is nested in capture batch here. No disease effects computed."}
    Path("reports/discovery_design_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = ["# Discovery metadata-only design audit", "", "These are candidate designs, not frozen models. No expression values or disease effects were read.", "", "| Candidate set | Terms in addition to intercept and disease | n | Rank / columns | Residual df | Structural gate |", "|---|---|---:|---:|---:|---|"]
    for name, items in results.items():
        for r in items:
            lines.append(f"| {name} | {', '.join(r['terms'])} | {r['n']} | {r['rank']} / {r['n_columns']} | {r['residual_df']} | {r['structural_gate_pass']} |")
    lines += ["", "The author-adjudicated phase scenario contains12 ENDO and10 non-ENDO donors. The strict histology-concordance scenario contains11 and9. All22 meet>=50 cells under the externally derived fibroblast definition, which remains subject to our annotation crosswalk. Under author-adjudicated phase, both disease groups coexist in capture batches12,14,16 and18; under strict histology concordance they coexist only in14 and16. The remaining batches contribute nuisance-parameter/variance information; they do not create within-batch disease comparisons. A full-rank design is not evidence that confounding is absent.", "", "Do not drop capture-batch adjustment simply because the smaller model gives more residual degrees of freedom. Sequencing batch is nested in capture batch; adding both produces redundant columns. Adenomyosis-adjusted and adenomyosis-free analyses answer different questions. Under author-adjudicated phase, excluding adenomyosis leaves8 cases and10 controls; structural results are shown above. Under strict histology concordance, it leaves only7 cases and9 controls and fails the fixed8-per-group minimum, even if a numerical regression could run. When adenomyosis is constant after exclusion, its redundant model column must be omitted, with this rule stated before effects.", "", "Machine-readable inputs, cohort IDs, covariate columns, source hashes, structural checks and batch tables are in discovery_design_audit.json. The script checks a known perfectly confounded artificial metadata copy as a negative control; it does not simulate biological results.", "", "YELLOW: understanding not yet assessed. Explain why an identifiable adjusted disease coefficient can still rely on few mixed batches, and why this is neither external validation nor a guarantee of power."]
    Path("reports/discovery_design_audit.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({k: [{a:r[a] for a in ['terms','n','rank','n_columns','residual_df','structural_gate_pass']} for r in v] for k,v in results.items()}, indent=2))

if __name__ == "__main__":
    main()
