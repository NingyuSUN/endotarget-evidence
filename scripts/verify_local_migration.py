#!/usr/bin/env python3
"""Verify migrated historical EndoTarget evidence without fitting or changing it.

Only /home/analyst/te-prefixed historical paths are mapped to --root for reads.
Current delivery maps select v1 outputs and v2 inputs intentionally: old v1
input hashes reference historical plotting code, not the current implementation.
"""
import argparse
import csv
import hashlib
import json
import math
import pathlib
import re
import subprocess
import sys
from datetime import datetime, timezone

from build_replication_report import (
    bootstrap_summary, lodo_summary, rows, truth, validate_counts,
)

SOURCE_COMMIT = "4cc7580851cbfcb776385ac0ae858fe5f81e818b"
HISTORICAL_ROOT = pathlib.Path("/home/analyst/te")
FREEZES = {
    "receipts/statistical_pre_effect_freeze.json":
        ("STATISTICAL_IMPLEMENTATION_FROZEN_PRE_EFFECT", 43),
    "receipts/discovery_rank_freeze.json":
        ("DISCOVERY_FROZEN_PRE_EXTERNAL", 74),
    "receipts/replication_implementation_freeze.json":
        ("REPLICATION_IMPLEMENTATION_FROZEN_PRE_EXTERNAL", 95),
}
METHODS = ("B0", "B1", "G1")


def require(condition, message):
    if not condition:
        raise ValueError(message)


class IntegrityChecker:
    def __init__(self, root):
        self.root = pathlib.Path(root).resolve(strict=True)
        self.cache = {}
        self.hash_checks = 0
        self.historical_paths = set()
        self.map_checks = {}

    @property
    def unique_files_hashed(self):
        return len(self.cache)

    def resolve(self, name, base=None):
        require(isinstance(name, (str, pathlib.Path)), "Invalid artifact path")
        name = str(name)
        require(bool(name) and "\x00" not in name, "Invalid artifact path")
        p = pathlib.Path(name)
        boundary = self.root if base is None else self.resolve(base)
        if p.is_absolute():
            if p.is_relative_to(HISTORICAL_ROOT):
                self.historical_paths.add(name)
                p = self.root / p.relative_to(HISTORICAL_ROOT)
            elif not p.is_relative_to(self.root):
                raise ValueError("Unsupported absolute path: " + name)
        else:
            p = boundary / p
        p = p.resolve()
        require(p.is_relative_to(self.root), "Artifact escapes local project: " + name)
        require(p.is_relative_to(boundary), "Artifact escapes selected directory: " + name)
        return p

    def digest(self, name, base=None):
        p = self.resolve(name, base)
        require(p.is_file(), "Missing artifact: " + str(p))
        stat = p.stat()
        stamp = (stat.st_size, stat.st_mtime_ns, stat.st_ino)
        if p in self.cache:
            require(self.cache[p][1] == stamp, "Artifact changed during verification: " + str(p))
            return self.cache[p][0]
        digest = hashlib.sha256()
        with p.open("rb") as handle:
            for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
                digest.update(block)
        after = p.stat()
        require(stamp == (after.st_size, after.st_mtime_ns, after.st_ino),
                "Artifact changed during hashing: " + str(p))
        self.cache[p] = (digest.hexdigest(), stamp, stat.st_size)
        return digest.hexdigest()

    def check_hash(self, name, expected, base=None):
        if isinstance(expected, dict):
            expected = expected.get("sha256")
        require(isinstance(expected, str) and re.fullmatch("[0-9a-f]{64}", expected),
                "Invalid expected SHA256: " + str(name))
        actual = self.digest(name, base)
        self.hash_checks += 1
        require(actual == expected, "SHA256 mismatch: " + str(name))

    def load(self, name):
        p = self.resolve(name)
        require(p.is_file(), "Missing artifact: " + str(p))
        return json.loads(p.read_text())

    def check_map(self, label, mapping, base=None):
        require(isinstance(mapping, dict) and bool(mapping), "Empty hash map: " + label)
        for name, expected in mapping.items():
            self.check_hash(name, expected, base)
        self.map_checks[label] = len(mapping)
        return len(mapping)


def git(root, *arguments):
    return subprocess.check_output(["git", "-C", str(root), *arguments], text=True).strip()


def verify_historical_receipt(checker, name):
    """Anchor existing receipt bytes to the selected preserved source commit."""
    p = checker.resolve(name)
    relative = str(p.relative_to(checker.root))
    original = subprocess.check_output(
        ["git", "-C", str(checker.root), "show", SOURCE_COMMIT + ":" + relative]
    )
    checker.check_hash(name, hashlib.sha256(original).hexdigest())


def verify_freezes(checker):
    result = []
    union = set()
    for filename, (phase, count) in FREEZES.items():
        verify_historical_receipt(checker, filename)
        obj = checker.load(filename)
        require(obj.get("phase") == phase, "Unexpected freeze phase: " + filename)
        require(len(obj["artifacts"]) == count, "Changed freeze coverage: " + filename)
        checker.check_map(filename + ":artifacts", obj["artifacts"])
        union.update(str(checker.resolve(p).relative_to(checker.root)) for p in obj["artifacts"])
        result.append({"receipt": filename, "phase": phase, "artifacts": count})
    require(len(union) == 95, "Changed frozen artifact union")
    return {"manifests": result, "unique_artifacts": len(union)}


def verify_delivery(checker):
    name = "receipts/scientific_delivery_verification.json"
    verify_historical_receipt(checker, name)
    delivery = checker.load(name)
    for specification, expected_count in delivery["hash_map_checks"].items():
        filename, field = specification.split(":", 1)
        verify_historical_receipt(checker, filename)
        obj = checker.load(filename)
        require(len(obj[field]) == expected_count, "Changed map coverage: " + specification)
        checker.check_map(specification, obj[field])
    checker.check_map(name + ":delivery_sha256", delivery["delivery_sha256"])
    final_name = "receipts/final_independent_delivery_review.json"
    verify_historical_receipt(checker, final_name)
    review = checker.load(final_name)
    checker.check_map(final_name + ":evidence_sha256", review["evidence_sha256"])
    # Primary execution receipts are preserved even where their results are inline.
    verify_historical_receipt(checker, "receipts/discovery_primary_execution.json")
    return {
        "delivery_map_specifications": len(delivery["hash_map_checks"]),
        "final_review_evidence_entries": len(review["evidence_sha256"]),
        "historical_v1_policy": "verify v1 outputs and current v2 inputs as selected by original delivery map",
    }


def verify_bootstrap(checker):
    directory = "work/statistical/bootstrap"
    summary, table = bootstrap_summary(checker.resolve(directory), False)
    require(summary is not None, "Missing bootstrap summary")
    checker.check_hash(directory + "/binding.rds", summary["binding_sha256"])
    frozen = checker.load(directory + "/draw_freeze.json")
    require(frozen["ledger_saved_before_any_fit"] is True, "Bootstrap ledger timing changed")
    checker.check_hash(directory + "/draw_ledger.tsv", frozen["ledger_sha256"])
    checker.check_map(directory + "/draw_freeze.json:inputs", frozen["binding"]["inputs"])
    mask = rows(checker.resolve("work/statistical/heca_primary/evaluation/bootstrap_interface.tsv"))
    family = len(frozen["binding"]["family"])
    expected_rates = {
        method + "_top" + str(k): sum(
            r["method"] == method and int(r["k"]) == k and truth(r["original_evaluable"])
            for r in mask
        ) for method in METHODS for k in (20, 50)
    }
    actual_files = {p.name for p in checker.resolve(directory).glob("replicate_*.json")}
    expected_files = {f"replicate_{i:03d}.json" for i in range(1, 501)}
    require(actual_files == expected_files, "Bootstrap draw file coverage changed")
    actual_sidecars = {p.name for p in checker.resolve(directory).glob("replicate_*.json.sha256")}
    require(actual_sidecars == {p + ".sha256" for p in expected_files},
            "Bootstrap sidecar coverage changed")
    valid_count = 0
    table_by_draw = {int(row["replicate"]): row for row in table}
    for number in range(1, 501):
        filename = directory + f"/replicate_{number:03d}.json"
        sidecar = checker.resolve(filename + ".sha256")
        require(sidecar.is_file(), "Missing bootstrap checksum: " + str(sidecar))
        checksum = sidecar.read_text().splitlines()
        require(len(checksum) == 1, "Invalid bootstrap checksum sidecar")
        checker.check_hash(filename, checksum[0])
        draw = checker.load(filename)
        require(draw["replicate"] == number and type(draw["valid"]) is bool,
                "Unsafe bootstrap draw identity")
        require(draw["binding_sha256"] == summary["binding_sha256"],
                "Bootstrap draw binding mismatch")
        require(draw["paired_replicate_valid"] is draw["valid"],
                "Bootstrap paired validity mismatch")
        row = table_by_draw[number]
        require(truth(row["valid"]) is draw["valid"], "Bootstrap table/draw validity mismatch")
        require(row["reason"] == draw["reason"], "Bootstrap table/draw reason mismatch")
        if draw["valid"]:
            valid_count += 1
            require(draw["family_size"] == family, "Bootstrap fixed family changed")
            require(set(draw["rates"]) == set(expected_rates), "Bootstrap rate keys changed")
            for key, denominator in expected_rates.items():
                rate = draw["rates"][key]
                k = int(key.split("_top")[1])
                agree = rate["agree"]
                require(type(agree) in (int, float) and math.isfinite(agree)
                        and agree == int(agree) and 0 <= agree <= denominator,
                        "Invalid bootstrap agreement count")
                require(rate["selected"] == k and rate["evaluable"] == denominator,
                        "Bootstrap fixed denominator changed")
                require(denominator > 0 and math.isfinite(rate["rate"])
                        and abs(rate["rate"] - agree / denominator) < 1e-14,
                        "Bootstrap fixed rate changed")
            difference = draw["rates"]["G1_top20"]["rate"] - draw["rates"]["B0_top20"]["rate"]
            require(math.isfinite(draw["difference"])
                    and abs(draw["difference"] - difference) < 1e-14,
                    "Bootstrap paired difference changed")
            require(abs(float(row["G1_minus_B0_top20"]) - difference) < 1e-14,
                    "Bootstrap table/draw difference mismatch")
        else:
            require(isinstance(draw["reason"], str) and bool(draw["reason"]),
                    "Failed bootstrap draw lacks reason")
    require(valid_count == summary["valid"], "Bootstrap valid draw total changed")
    return {"draws": 500, "sidecar_SHA256_checks": 500,
            "valid": valid_count, "invalid_preserved": 500 - valid_count,
            "conditional_interval": summary["interval"], "fits_executed": 0}


def verify_lodo(checker):
    directory = "work/statistical/discovery_lodo"
    summary, metrics = lodo_summary(checker.resolve(directory), False)
    require(summary is not None, "Missing LODO summary")
    checker.check_hash(directory + "/binding.json", summary["binding_sha256"])
    binding = checker.load(directory + "/binding.json")
    checker.check_map(directory + "/binding.json:sha256", binding["sha256"])
    donors = rows(checker.resolve("manifests/final_discovery_clinical_selected.tsv"))
    require(len(donors) == 22, "Changed discovery donor count")
    require({p.name for p in checker.resolve(directory).glob("fold_*") if p.is_dir()}
            == {f"fold_{i:02d}" for i in range(1, 23)}, "LODO fold coverage changed")
    row_by_key = {(int(row["fold"]), row["method"]): row for row in metrics}
    output_checks = 0
    for number, donor in enumerate(donors, 1):
        fold = directory + f"/fold_{number:02d}"
        checkpoint = checker.load(fold + "/checkpoint.json")
        require(checkpoint["binding_sha256"] == summary["binding_sha256"]
                and checkpoint["dropped_donor"] == donor["donor_id"], "LODO checkpoint identity changed")
        require(checkpoint["valid"] is True and not checkpoint["failed_fold_replaced"]
                and checkpoint["fit_exit_status"] == 0, "LODO completion record changed")
        output_checks += checker.check_map(
            fold + "/checkpoint.json:output_sha256", checkpoint["output_sha256"], base=fold
        )
        actual = {
            str(p.relative_to(checker.resolve(fold)))
            for p in checker.resolve(fold).rglob("*")
            if p.is_file() and p.name != "checkpoint.json"
        }
        require(actual == set(checkpoint["output_sha256"]), "LODO checkpoint output coverage changed")
        for method in METHODS:
            saved = checkpoint["metrics"][method]
            row = row_by_key[(number, method)]
            require(truth(row["valid"]) and row["dropped_donor"] == donor["donor_id"],
                    "LODO metrics donor/validity changed")
            require(int(row["intersection"]) == saved["intersection"]
                    and abs(float(row["retention"]) - saved["retention"]) < 1e-14
                    and abs(float(row["jaccard"]) - saved["jaccard"]) < 1e-14,
                    "LODO checkpoint/summary metrics mismatch")
    return {"folds": 22, "valid_folds": summary["valid_folds"],
            "fold_output_SHA256_checks": output_checks, "metrics_rows": len(metrics),
            "fits_executed": 0}


def verify_reported_results(checker):
    evaluation = checker.load("work/statistical/heca_primary/evaluation/summary.json")
    selected = rows(checker.resolve("work/statistical/heca_primary/evaluation/selected_genes.tsv"))
    validate_counts(evaluation, selected)
    observed = {
        method: {
            str(k): [evaluation["methods"][method][str(k)][key]
                     for key in ("agree", "evaluable", "selected")]
            for k in (20, 50)
        } for method in METHODS
    }
    expected = {
        "B0": {"20": [14, 20, 20], "50": [36, 49, 50]},
        "B1": {"20": [3, 13, 20], "50": [16, 35, 50]},
        "G1": {"20": [6, 19, 20], "50": [26, 47, 50]},
    }
    require(observed == expected, "Historical concordance changed")
    tested = {}
    for cohort, count in (("discovery", 17733), ("heca", 17643)):
        table = rows(checker.resolve(f"work/statistical/{cohort}_primary/de/de_all.tsv"))
        require(bool(table) and "adj.P.Val" in table[0], "Missing DE adjusted-P column")
        significant = sum(float(row["adj.P.Val"]) < .05 for row in table
                          if row["adj.P.Val"] not in ("", "NA"))
        require(len(table) == count and significant == 0, "Historical DE count/significance changed")
        tested[cohort] = {"tested_genes": len(table), "BH_q_lt_05": significant}
    discovery = checker.load("receipts/discovery_primary_execution.json")
    require(discovery["expression_pass"] == 17733 and discovery["BH_q_lt_05"] == 0,
            "Discovery execution receipt count changed")
    report = checker.load("receipts/replication_v2_report.json")
    require(report["bootstrap"]["valid"] == 487 and report["bootstrap"]["invalid"] == 13
            and report["lodo"]["valid_folds"] == 22, "Report resampling counts changed")
    return {"concordance_agree_evaluable_selected": observed, "DE": tested,
            "selected_evaluation_rows": len(selected), "scientific_results_recomputed": False}


def run(root):
    checker = IntegrityChecker(root)
    report = {
        "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_commit": SOURCE_COMMIT, "historical_root": str(HISTORICAL_ROOT),
        "local_root": str(checker.root), "local_branch": git(checker.root, "branch", "--show-current"),
        "local_head_before_receipt": git(checker.root, "rev-parse", "HEAD"),
        "scientific_validation": False, "therapeutic_validation": False,
        "container_validation": False, "scientific_fits_executed": 0,
        "scope": "frozen and current delivery evidence plus all bootstrap and LODO checkpoints; not full raw import or a whole-transfer checksum comparison",
        "checks": {}, "errors": [],
    }
    subprocess.check_call(["git", "-C", str(checker.root), "merge-base", "--is-ancestor",
                           SOURCE_COMMIT, report["local_head_before_receipt"]])
    stages = [("freezes", verify_freezes), ("delivery", verify_delivery),
              ("bootstrap", verify_bootstrap), ("LODO", verify_lodo),
              ("reported_results", verify_reported_results)]
    for label, function in stages:
        try:
            report["checks"][label] = {"status": "PASS", **function(checker)}
        except (ValueError, OSError, KeyError, TypeError, subprocess.CalledProcessError) as error:
            report["checks"][label] = {"status": "FAIL"}
            report["errors"].append({"stage": label, "error": str(error)})
    verifier_sha256 = checker.digest(pathlib.Path(__file__).resolve())
    report.update({
        "status": "PASS_MIGRATION_EVIDENCE_INTEGRITY" if not report["errors"] else "FAIL_MIGRATION_EVIDENCE_INTEGRITY",
        "hash_checks": checker.hash_checks,
        "unique_files_hashed": checker.unique_files_hashed,
        "unique_bytes_hashed": sum(entry[2] for entry in checker.cache.values()),
        "historical_absolute_paths_mapped": sorted(checker.historical_paths),
        "completed_hash_maps": checker.map_checks,
        "verifier_sha256": verifier_sha256,
    })
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    report = run(args.root)
    output = args.root.resolve() / "receipts/local_migration_integrity.json"
    temporary = output.with_suffix(".json.partial")
    temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    temporary.replace(output)
    print(json.dumps({key: report[key] for key in
                     ("status", "unique_files_hashed", "unique_bytes_hashed", "hash_checks", "errors")}))
    return 0 if not report["errors"] else 1


if __name__ == "__main__":
    sys.exit(main())
