#!/usr/bin/env python3
"""Discovery-only B0/B1/G1 on one induced graph; no external data interface."""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
import math
import pathlib
import platform
import re
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
GENE = re.compile(r"^ENSG[0-9]{11}$")
DE_FIELDS = ["gene_id", "logFC", "t", "P.Value", "adj.P.Val", "AveExpr"]
RANK_FIELDS = ["rank", "gene_id", "method", "score", "logFC", "discovery_sign", "t", "P.Value", "adj.P.Val", "AveExpr", "weighted_degree"]
REFERENCE_SHA256 = "84a4313292da6d6d69b76a4f721893aa89775beebeed9eabcb263fe93e797033"
ALGORITHM = {
    "schema_version": 1,
    "minimum_universe_genes": 50,
    "universe": "complete_discovery_expression_pass_de_all_induced_graph_then_remove_zero_degree_and_require_all_U_fits",
    "stable_gene_id_pattern": "^ENSG[0-9]{11}$",
    "reference_sha256": REFERENCE_SHA256,
    "reference_kind": "STRING12.5_direct_experiments_database_functional_association",
    "reference_schema": ["gene1", "gene2", "weight"],
    "reference_weight": "unchanged_max_direct_score_divided_by_1000",
    "minimum_edge_weight": 0.7,
    "maximum_edge_weight": 1.0,
    "duplicate_reference_edges": "reject_already_collapsed_immutable_reference",
    "self_loops": "reject_already_removed_immutable_reference",
    "B0": "descending_abs_moderated_discovery_t",
    "B1": "descending_weighted_degree_same_U",
    "G1": "restart_random_walk_symmetric_adjacency_column_normalized",
    "seed": "abs_discovery_t_normalized_to_sum_1_on_U",
    "restart_probability": 0.5,
    "l1_tolerance": 1e-10,
    "max_iterations": 10000,
    "tie_break": "ascending_stable_gene_id_exact_score_ties",
    "direction": "sign_discovery_logFC_independently_preserved_zero_is_zero",
    "top_k": [20, 50],
    "external_data_used": False,
    "parameter_search": False,
}


def digest(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def table_rows(path, fields, exact=False):
    """Read TSV/gzip TSV, rejecting malformed headers and ragged rows."""
    path = pathlib.Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames
        if not header or len(header) != len(set(header)) or not set(fields).issubset(header):
            raise ValueError(f"Invalid or duplicate TSV schema in {path}")
        if exact and header != fields:
            raise ValueError(f"Reference schema must be exactly {fields}")
        for line, item in enumerate(reader, 2):
            if None in item or any(v is None for v in item.values()):
                raise ValueError(f"Ragged TSV row {line} in {path}")
            yield item


def validated_de(rows):
    result = {}
    for i, source in enumerate(rows, 1):
        if not all(f in source for f in DE_FIELDS):
            raise ValueError(f"Discovery DE row {i} is missing required schema")
        gid = source["gene_id"]
        if not isinstance(gid, str) or not GENE.fullmatch(gid):
            raise ValueError(f"Invalid stable gene_id in discovery row {i}")
        if gid in result:
            raise ValueError(f"Duplicate discovery gene_id {gid}")
        item = {"gene_id": gid}
        for field in DE_FIELDS[1:]:
            try:
                number = float(source[field])
            except (ValueError, TypeError, OverflowError) as exc:
                raise ValueError(f"Invalid discovery numeric field {gid}/{field}") from exc
            if not math.isfinite(number):
                raise ValueError(f"Nonfinite discovery numeric field {gid}/{field}")
            if field in ("P.Value", "adj.P.Val") and not 0 <= number <= 1:
                raise ValueError(f"Discovery probability outside [0,1]: {gid}/{field}")
            item[field] = number
        result[gid] = item
    if not result:
        raise ValueError("Discovery DE has no fitted rows")
    return result


def read_de(path):
    return list(validated_de(table_rows(path, DE_FIELDS)).values())


def induced_edges(rows, expressed):
    """Validate whole immutable reference before retaining discovery-induced edges."""
    seen = set()
    kept = []
    count = 0
    for row in rows:
        if set(row) != {"gene1", "gene2", "weight"}:
            raise ValueError("Invalid reference edge schema")
        a, b = row["gene1"], row["gene2"]
        if not isinstance(a, str) or not isinstance(b, str) or not GENE.fullmatch(a) or not GENE.fullmatch(b):
            raise ValueError("Invalid reference stable gene identifier")
        if a == b:
            raise ValueError("Reference self loop in immutable collapsed graph")
        pair = tuple(sorted((a, b)))
        if pair in seen:
            raise ValueError("Duplicate undirected reference edge; immutable graph must already be max-collapsed")
        seen.add(pair)
        try:
            weight = float(row["weight"])
        except (ValueError, TypeError, OverflowError) as exc:
            raise ValueError("Invalid reference edge weight") from exc
        if not math.isfinite(weight) or not .7 <= weight <= 1:
            raise ValueError("Reference weight must be finite and within [0.7,1]")
        count += 1
        if a in expressed and b in expressed:
            kept.append((pair[0], pair[1], weight))
    return sorted(kept), count


def propagate(seed, left, right, weights, degree, *, max_iterations=10000):
    """Sparse W*p by two scatter sums; O(nodes + edges) memory, no dense matrix."""
    seed = np.asarray(seed, dtype=np.float64)
    degree = np.asarray(degree, dtype=np.float64)
    left, right = np.asarray(left, dtype=np.int64), np.asarray(right, dtype=np.int64)
    weights = np.asarray(weights, dtype=np.float64)
    if seed.ndim != 1 or len(seed) == 0 or degree.shape != seed.shape:
        raise ValueError("Invalid propagation node arrays")
    if left.ndim != 1 or left.shape != right.shape or left.shape != weights.shape:
        raise ValueError("Invalid propagation edge arrays")
    if not isinstance(max_iterations, int) or isinstance(max_iterations, bool) or max_iterations < 1:
        raise ValueError("Invalid propagation iteration count")
    if not np.all(np.isfinite(seed)) or not np.all(seed >= 0) or not math.isclose(float(seed.sum()), 1., abs_tol=1e-12):
        raise ValueError("Propagation seed must be finite nonnegative and sum to 1")
    if not np.all(np.isfinite(degree)) or not np.all(degree > 0):
        raise ValueError("Propagation degree must be finite and positive")
    if not np.all(np.isfinite(weights)) or not np.all(weights > 0):
        raise ValueError("Propagation weights must be finite and positive")
    if np.any(left < 0) or np.any(right < 0) or np.any(left >= len(seed)) or np.any(right >= len(seed)):
        raise ValueError("Propagation edge index outside nodes")
    observed = np.bincount(left, weights=weights, minlength=len(seed)) + np.bincount(right, weights=weights, minlength=len(seed))
    if not np.allclose(observed, degree, rtol=1e-14, atol=1e-14):
        raise ValueError("Propagation degree does not match symmetric graph")
    n = len(seed)
    def walk(p):
        # Column j transfers p[j] to each neighbor i as A[i,j] / degree[j].
        return (np.bincount(left, weights=weights * p[right] / degree[right], minlength=n)
                + np.bincount(right, weights=weights * p[left] / degree[left], minlength=n))
    p = seed.copy()
    maximum_mass_error = abs(float(p.sum()) - 1.)
    for iteration in range(1, max_iterations + 1):
        following = .5 * seed + .5 * walk(p)
        error = abs(float(following.sum()) - 1.)
        maximum_mass_error = max(maximum_mass_error, error)
        if not np.all(np.isfinite(following)) or np.any(following < 0) or error > 1e-12:
            raise ValueError("Propagation probability mass/nonnegative invariant failed")
        change = float(np.abs(following - p).sum())
        p = following
        if change < 1e-10:
            residual = float(np.abs(.5 * seed + .5 * walk(p) - p).sum())
            return p, {"iterations": iteration, "final_l1_change": change,
                       "fixed_point_l1_residual": residual, "mass_sum": float(p.sum()),
                       "maximum_mass_error": maximum_mass_error, "converged": True}
    raise ValueError(f"Propagation failed to converge after {max_iterations} iterations")


def rank(de_rows, edge_rows, *, eligibility_rows=None, test_min_universe=None):
    """Pure input-data API. A named test override exists only for tiny hand-solvable graphs."""
    minimum = 50 if test_min_universe is None else test_min_universe
    if not isinstance(minimum, int) or isinstance(minimum, bool) or minimum < 1:
        raise ValueError("Invalid test minimum universe count")
    de = validated_de(de_rows)
    expressed = set()
    for row in (de.values() if eligibility_rows is None else eligibility_rows):
        gid = row.get("gene_id")
        if not isinstance(gid, str) or not GENE.fullmatch(gid):
            raise ValueError("Invalid complete discovery expression-pass gene identifier")
        if gid in expressed:
            raise ValueError("Duplicate complete discovery expression-pass gene identifier")
        expressed.add(gid)
    if not expressed:
        raise ValueError("Complete discovery expression-pass family is empty")
    if set(de) - expressed:
        raise ValueError("Fitted discovery genes outside complete expression-pass family")
    edges, reference_count = induced_edges(edge_rows, expressed)
    universe = sorted({g for a, b, _ in edges for g in (a, b)})
    missing = sorted(set(universe) - set(de))
    if missing:
        raise ValueError(f"Induced-U genes missing finite discovery fit: {len(missing)}; examples {missing[:5]}")
    if len(universe) < minimum:
        raise ValueError(f"Final induced nonisolated universe has {len(universe)} genes; requires >= {minimum}")
    index = {g: i for i, g in enumerate(universe)}
    left = np.array([index[a] for a, b, w in edges], dtype=np.int64)
    right = np.array([index[b] for a, b, w in edges], dtype=np.int64)
    weights = np.array([w for a, b, w in edges], dtype=np.float64)
    degree = np.bincount(left, weights=weights, minlength=len(universe)) + np.bincount(right, weights=weights, minlength=len(universe))
    abs_t = np.array([abs(de[g]["t"]) for g in universe], dtype=np.float64)
    largest = float(abs_t.max())
    if largest == 0:
        raise ValueError("All discovery seeds on final universe are zero")
    # Scaling first avoids overflow when valid finite t statistics sum above float max.
    scaled = abs_t / largest
    seed = scaled / scaled.sum()
    propagated, diagnostic = propagate(seed, left, right, weights, degree)
    rankings = {}
    for method, scores in [("B0", abs_t), ("B1", degree), ("G1", propagated)]:
        order = sorted(range(len(universe)), key=lambda i: (-float(scores[i]), universe[i]))
        table = []
        for position, i in enumerate(order, 1):
            source = de[universe[i]]
            table.append({"rank": position, **source, "method": method, "score": float(scores[i]),
                          "discovery_sign": (source["logFC"] > 0) - (source["logFC"] < 0),
                          "weighted_degree": float(degree[i])})
        rankings[method] = table
    return {"universe": universe, "edges": edges, "rankings": rankings,
            "excluded": [{"gene_id": g, "reason": "zero_degree_in_discovery_induced_graph"} for g in sorted(expressed - set(universe))],
            "counts": {"discovery_expressed_genes": len(expressed), "finite_discovery_fit_genes": len(de), "reference_edges": reference_count,
                       "induced_edges": len(edges), "final_universe_genes": len(universe),
                       "excluded_induced_zero_degree": len(expressed) - len(universe)},
            "propagation": diagnostic, "test_min_universe_override": test_min_universe}


def write_table(path, rows, fields):
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({f: format(row[f], ".17g") if isinstance(row[f], float) else row[f] for f in fields})


def execute(de_path, edge_path, output_dir, config_path, *, de_all_path=None, synthetic_test=False):
    config_path, de_path, edge_path, output_dir = map(pathlib.Path, (config_path, de_path, edge_path, output_dir))
    if de_all_path is None and not synthetic_test:
        raise ValueError("Production ranking requires --de-all complete expression-pass family")
    de_all_path = pathlib.Path(de_path if de_all_path is None else de_all_path)
    config = json.loads(config_path.read_text())
    if config != ALGORITHM:
        raise ValueError("Ranking algorithm configuration differs from frozen implementation")
    hashes = {"discovery_de": digest(de_path), "discovery_de_all": digest(de_all_path), "reference_edges": digest(edge_path), "configuration": digest(config_path), "code": digest(__file__)}
    if not synthetic_test and hashes["reference_edges"] != config["reference_sha256"]:
        raise ValueError("Pinned reference edge hash mismatch")
    partial = output_dir.with_name(output_dir.name + ".partial")
    if output_dir.exists() or partial.exists():
        raise ValueError("Existing ranking output/partial preserved; use a new output directory")
    result = rank(read_de(de_path), table_rows(edge_path, ["gene1", "gene2", "weight"], exact=True), eligibility_rows=table_rows(de_all_path, DE_FIELDS))
    for name, path in [("discovery_de", de_path), ("discovery_de_all", de_all_path), ("reference_edges", edge_path), ("configuration", config_path), ("code", pathlib.Path(__file__))]:
        if digest(path) != hashes[name]:
            raise ValueError(f"Input changed during ranking: {name}")
    partial.mkdir(parents=True, exist_ok=False)
    for method, table in result["rankings"].items():
        write_table(partial / f"{method}.tsv", table, RANK_FIELDS)
        for k in (20, 50):
            write_table(partial / f"{method}_top{k}.tsv", table[:k], RANK_FIELDS)
    write_table(partial / "universe.tsv", [{"gene_id": g} for g in result["universe"]], ["gene_id"])
    write_table(partial / "excluded_genes.tsv", result["excluded"], ["gene_id", "reason"])
    write_table(partial / "induced_edges.tsv", [{"gene1": a, "gene2": b, "weight": w} for a, b, w in result["edges"]], ["gene1", "gene2", "weight"])
    outputs = {p.name: {"sha256": digest(p), "bytes": p.stat().st_size} for p in sorted(partial.iterdir())}
    provenance = {
        "status": "SYNTHETIC_TEST_ONLY_NOT_A_SCIENTIFIC_FREEZE" if synthetic_test else "DISCOVERY_RANKING_ARTIFACT_NOT_COMPLETE_EXTERNAL_FREEZE",
        "finished_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "inputs": {name: {"path": str(path.resolve()), "sha256": hashes[name]} for name, path in [("discovery_de", de_path), ("discovery_de_all", de_all_path), ("reference_edges", edge_path), ("configuration", config_path)]},
        "implementation": {"path": str(pathlib.Path(__file__).resolve()), "sha256": hashes["code"]},
        "software": {"python": platform.python_version(), "numpy": np.__version__, "executable": sys.executable},
        "algorithm": config, "counts": result["counts"], "propagation": result["propagation"],
        "outputs": outputs, "universe_sha256": outputs["universe.tsv"]["sha256"],
        "induced_graph_sha256": outputs["induced_edges.tsv"]["sha256"],
        "external_data_used": False, "test_min_universe_override": None,
        "scientific_validation": False, "complete_external_freeze": False,
        "direction_meaning": "discovery_expression_sign_not_intervention_activation_or_inhibition",
        "reference_meaning": "functional_association_not_physical_binding_or_causal_regulation",
    }
    (partial / "provenance.json").write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
    partial.rename(output_dir)
    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--de", required=True, type=pathlib.Path, help="Finite discovery fitted rows")
    parser.add_argument("--de-all", type=pathlib.Path, help="Complete expression-pass DE table, required for production; failed fitted fields may be NA")
    parser.add_argument("--edges", type=pathlib.Path, default=ROOT / "data/reference/string/direct_evidence_gene_edges.tsv.gz")
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    parser.add_argument("--config", type=pathlib.Path, default=ROOT / "config/ranking_algorithm.json")
    parser.add_argument("--synthetic-test", action="store_true", help="Permit a synthetic unpinned graph; keeps >=50 gate and marks all provenance synthetic")
    args = parser.parse_args()
    try:
        provenance = execute(args.de, args.edges, args.output_dir, args.config, de_all_path=args.de_all, synthetic_test=args.synthetic_test)
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        parser.exit(1, f"Ranking failed: {exc}\n")
    print(json.dumps({"status": provenance["status"], "output_dir": str(args.output_dir), "counts": provenance["counts"], "propagation": provenance["propagation"]}, indent=2))


if __name__ == "__main__":
    main()
