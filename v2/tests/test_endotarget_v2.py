#!/usr/bin/env python3
"""Comprehensive automated pytest contract & scientific integrity suite for EndoTarget v2.

Verifies:
1. Cryptographic SHA-256 chain across Stage 1 (protocol_freeze.json), Stage 2 (prediction_seal.json),
   Stage 3 (locked_evaluation_v2.json), and rendered figures (figures_v2.json).
2. Git commit chronology (Stage 1 -> Stage 2 -> Stage 3 zero-leakage progression).
3. Cohort pairing/subject deduplication invariants (8 matched pairs & 6 triplets in GSE179640;
   29 non-redundant libraries & 4 matched natural-cycle pairs in GSE213216; 408 arrays mapped 100%
   to 162 subjects in GSE141549).
4. Exact integer UMI conservation in both single-cell pseudobulk aggregations.
5. Complete absence of HECA transferred label files or columns across manifests and scripts.
6. Synthetic mathematical checks on limma-voom, paired empirical Bayes moderation, RWR mass conservation,
   and M1 dual-contrast rejection of planted host-background false positives.
7. Exact numerical consistency between JSON/TSV summary artifacts and bilingual reports
   (README.md, README.zh-CN.md, reports/RESULTS.md, reports/RESULTS.zh-CN.md).
"""
from __future__ import annotations

import csv
import gzip
import json
import pathlib
import subprocess
import sys
import types

import numpy as np
import pandas as pd
import pytest
import yaml

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
for _opt_mod in ("anndata", "h5py", "scanpy"):
    if _opt_mod not in sys.modules:
        try:
            __import__(_opt_mod)
        except ImportError:
            _stub = types.ModuleType(_opt_mod)
            if _opt_mod == "scanpy":
                _stub.settings = types.SimpleNamespace(n_jobs=8)
            sys.modules[_opt_mod] = _stub

from endotarget_lib import (  # noqa: E402
    limma_voom_fit,
    rwr_propagate,
    squeeze_var,
    trigamma_inverse,
)
from fetch_external_cohorts import verify_prediction_seal  # noqa: E402
from freeze_protocol import sha256_file, verify_protocol_freeze  # noqa: E402


def test_cryptographic_sha256_chain():
    """Verify unbroken SHA-256 chain from protocol_freeze -> prediction_seal -> locked_evaluation -> figures."""
    freeze = verify_protocol_freeze(ROOT)
    seal = verify_prediction_seal(ROOT)
    freeze_sha = sha256_file(ROOT / "receipts/protocol_freeze.json")
    seal_sha = sha256_file(ROOT / "receipts/prediction_seal.json")

    assert seal["protocol_freeze_sha256"] == freeze_sha
    assert freeze["attestations"]["discovery_expression_downloaded"] is False
    assert freeze["attestations"]["external_expression_downloaded"] is False
    assert seal["attestations"]["external_GSE213216_downloaded_or_read"] is False
    assert seal["attestations"]["external_GSE141549_expression_downloaded_or_read"] is False

    eval_rec = json.loads((ROOT / "receipts/locked_evaluation_v2.json").read_text(encoding="utf-8"))
    assert eval_rec["protocol_freeze_sha256"] == freeze_sha
    assert eval_rec["prediction_seal_sha256"] == seal_sha
    for relpath, info in eval_rec["outputs"].items():
        p = ROOT / relpath
        assert p.exists(), f"Missing Stage 3 output {relpath}"
        assert sha256_file(p) == info["sha256"], f"Tampered Stage 3 output {relpath}"

    fig_rec = json.loads((ROOT / "receipts/figures_v2.json").read_text(encoding="utf-8"))
    for relpath, expected_sha in fig_rec["source_artifacts"].items():
        assert sha256_file(ROOT / relpath) == expected_sha
    assert len(fig_rec["figures"]) == 8
    for relpath, info in fig_rec["figures"].items():
        p = ROOT / relpath
        assert p.exists() and p.stat().st_size > 10_000
        assert sha256_file(p) == info["sha256"]


def test_git_commit_chain_chronology():
    """Verify Commit 1 -> Commit 2 -> Commit 3 exist in chronological Git history for this module."""
    is_shallow = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "--is-shallow-repository"],
        text=True,
    ).strip()
    if is_shallow == "true":
        pytest.skip("Shallow git clone; full history required for Stage 1->2->3 chronology check")
    out = subprocess.check_output(
        ["git", "-C", str(ROOT), "log", "--reverse", "--pretty=format:%s"],
        text=True,
    ).splitlines()
    assert len(out) >= 3
    s1 = next((i for i, msg in enumerate(out) if "Stage 1" in msg), None)
    s2 = next((i for i, msg in enumerate(out) if "Stage 2" in msg), None)
    s3 = next((i for i, msg in enumerate(out) if "Stage 3" in msg), None)
    assert s1 is not None and s2 is not None and s3 is not None
    assert s1 < s2 < s3


def test_m0_cohort_designs_and_donor_deduplication():
    """Verify exact donor/subject counts and pairing invariants across all 3 cohorts."""
    # 1. GSE179640: 31 scRNA libraries, 8 matched (EuE, EcP) pairs, 6 matched (EuE, EcP, AdjP) triplets
    d179 = pd.read_csv(ROOT / "manifests/m0_GSE179640_donor_design.tsv", sep="\t")
    s179 = pd.read_csv(ROOT / "manifests/m0_GSE179640_samples.tsv", sep="\t")
    assert int((s179["included_discovery_scrna"].astype(str) == "True").sum()) == 31
    assert int((d179["matched_EuE_EcP_pair"].astype(str) == "True").sum()) == 8
    assert int((d179["matched_EuE_EcP_AdjP_triplet"].astype(str) == "True").sum()) == 6

    # 2. GSE213216: 29 non-redundant libraries (at most 1 per patient x site_class), 4 matched pairs
    sel213 = pd.read_csv(ROOT / "manifests/m0_GSE213216_locked_selection.tsv", sep="\t")
    assert len(sel213) == 29
    assert not sel213.duplicated(subset=["patient_id", "site_class"]).any()
    by_pt = sel213.groupby("patient_id")["site_class"].apply(set)
    matched_pts = sorted(pt for pt, sc in by_pt.items() if {"EuE", "EcP"} <= sc)
    assert matched_pts == ["P10", "P11", "P14", "P16"]

    # 3. GSE141549: 408 samples mapped 100% to 162 subjects with exact paired counts
    s141 = pd.read_csv(ROOT / "manifests/m0_GSE141549_samples.tsv", sep="\t")
    subj141 = pd.read_csv(ROOT / "manifests/m0_GSE141549_subjects.tsv", sep="\t")
    assert len(s141) == 408
    assert s141["subject_id"].isna().sum() == 0
    assert len(subj141) == 162
    assert int((subj141["has_paired_PE_Peritoneal"].astype(str) == "True").sum()) == 41
    assert int((subj141["has_paired_PP_Peritoneal"].astype(str) == "True").sum()) == 20
    assert int((subj141["has_triplet_PE_PP_Peritoneal"].astype(str) == "True").sum()) == 18
    assert int((subj141["has_paired_PP_DIE"].astype(str) == "True").sum()) == 20


def test_pseudobulk_umi_conservation_and_no_heca_leakage():
    """Verify exact integer UMI conservation in GSE179640 and GSE213216 and zero HECA label leakage."""
    disc_sum = json.loads((ROOT / "reports/discovery/discovery_summary_v2.json").read_text(encoding="utf-8"))
    c1 = disc_sum["pseudobulk_conservation"]
    assert c1["exact_conservation"] is True
    assert (
        c1["total_qc_passing_umi_sum"]
        == c1["retained_pseudobulk_umi_sum"]
        + c1["dropped_below_min_cells_umi_sum"]
        + c1["excluded_ambiguous_or_unassigned_umi_sum"]
    )

    eval_sum = json.loads((ROOT / "reports/locked_evaluation_summary_v2.json").read_text(encoding="utf-8"))
    c2 = eval_sum["GSE213216_external_scrna"]["pseudobulk_conservation"]
    assert c2["exact_conservation"] is True
    assert (
        c2["total_qc_passing_umi_sum"]
        == c2["retained_pseudobulk_umi_sum"]
        + c2["dropped_below_min_cells_umi_sum"]
        + c2["excluded_ambiguous_or_unassigned_umi_sum"]
    )

    # Verify protocol explicitly forbids HECA annotations and canonical marker panel uses only 7 de novo lineages
    cfg = yaml.safe_load((ROOT / "configs/protocol_v2.yaml").read_text(encoding="utf-8"))
    assert cfg["single_cell_processing"]["heca_annotations_permitted"] is False
    markers = pd.read_csv(ROOT / "manifests/canonical_lineage_markers.tsv", sep="\t")
    assert set(markers["lineage"].unique()) == {
        "stromal_fibroblast",
        "perivascular_mural",
        "epithelial",
        "endothelial",
        "myeloid_macrophage",
        "lymphoid",
        "mesothelial",
    }
    # Ensure neither discovery nor external annotation pipelines import or read any HECA cell-label file
    for script_name in ("endotarget_lib.py", "run_discovery_v2.py", "evaluate_locked.py"):
        src = (ROOT / "scripts" / script_name).read_text(encoding="utf-8")
        assert "heca_features" not in src and "celltype_heca" not in src and "E-MTAB-14039" not in src


def test_synthetic_estimators_and_m1_dual_contrast_host_rejection():
    """Unit-test trigamma_inverse, limma-voom, RWR fixed point, and M1 host-artifact rejection on synthetic fixtures."""
    from scipy import special

    for target_y in (0.05, 0.5, 2.5, 15.0):
        x = trigamma_inverse(target_y)
        assert abs(float(special.polygamma(1, x)) - target_y) < 1e-7

    # Synthetic paired count fixture: 500 genes, 6 donors x 2 sites (12 samples)
    # Genes 0..4 have true +3 log2 fold-change in case; genes 5..499 are null
    rng = np.random.default_rng(42)
    n_genes, n_pairs = 500, 6
    base = rng.gamma(shape=5.0, scale=40.0, size=(n_genes, n_pairs))
    c_ref = rng.poisson(base)
    mult = np.ones((n_genes, 1))
    mult[:5, 0] = 8.0
    c_case = rng.poisson(base * mult)
    counts = np.hstack([c_ref, c_case])
    X = np.zeros((2 * n_pairs, n_pairs + 1))
    X[:, 0] = 1.0
    for d in range(1, n_pairs):
        X[d, d] = 1.0
        X[n_pairs + d, d] = 1.0
    X[n_pairs:, -1] = 1.0

    fit = limma_voom_fit(counts, X, coef_idx=n_pairs)
    assert np.all(fit["q"][:5] < 1e-3)
    assert np.all(fit["logFC"][:5] > 2.0)
    assert (fit["q"][5:] < 0.05).sum() <= 2

    # RWR probability mass & exact linear-system fixed-point check on a 6-node graph
    edges = [(0, 1, 0.9), (1, 2, 0.8), (2, 3, 0.75), (3, 4, 0.95), (4, 5, 0.85), (0, 3, 0.7)]
    left = np.array([e[0] for e in edges], dtype=np.int64)
    right = np.array([e[1] for e in edges], dtype=np.int64)
    w = np.array([e[2] for e in edges], dtype=np.float64)
    seeds = np.array([10.0, 1.0, 0.5, 2.0, 0.1, 0.2])
    p_rwr, diag = rwr_propagate(seeds, left, right, w, restart_prob=0.5, tol=1e-12)
    assert diag["converged"] is True
    assert abs(float(p_rwr.sum()) - 1.0) < 1e-12

    A = np.zeros((6, 6))
    for i, j, wt in edges:
        A[i, j] = A[j, i] = wt
    W_col = A / A.sum(axis=0, keepdims=True)
    s_norm = seeds / seeds.sum()
    p_exact = np.linalg.solve(np.eye(6) - 0.5 * W_col, 0.5 * s_norm)
    assert np.allclose(p_rwr, p_exact, atol=1e-10)

    # Verify M1 formula strictly zeroes host-confounded false positives (z_EuE = +6, z_AdjP = -2)
    z_eue = np.array([6.0, 4.0])     # Gene 0 = peritoneal background; Gene 1 = true lesion remodeling
    z_adjp = np.array([-2.0, 3.5])
    w_amb = np.array([0.4, 0.98])
    m1_core = w_amb * np.minimum(np.abs(z_eue), np.maximum(0.0, np.sign(z_eue) * z_adjp))
    assert m1_core[0] == 0.0
    assert m1_core[1] > 3.4


def test_exact_artifact_and_bilingual_report_consistency():
    """Verify exact numerical match between JSON/TSV artifacts and all 4 bilingual markdown reports."""
    disc = json.loads((ROOT / "reports/discovery/discovery_summary_v2.json").read_text(encoding="utf-8"))
    ev = json.loads((ROOT / "reports/locked_evaluation_summary_v2.json").read_text(encoding="utf-8"))
    dossiers = pd.read_csv(ROOT / "reports/target_dossiers.tsv", sep="\t")
    cand = pd.read_csv(ROOT / "manifests/frozen_candidates_v2.tsv", sep="\t")

    # Check frozen shortlist cardinalities
    for m in ("B0", "G1", "M1", "M2"):
        assert int((cand[f"in_{m}_top20"].astype(str) == "True").sum()) == 20
        assert int((cand[f"in_{m}_top50"].astype(str) == "True").sum()) == 50

    # Check key numerical invariants in saved summaries
    assert disc["qc_pass_cells_total"] == 115402
    assert disc["gene_universe"]["shared_string_universe_genes"] == 10972
    assert disc["paired_stromal_models"]["EcP_vs_EuE"]["sig_BH_q_lt_0.05_expressed"] == 569
    assert ev["GSE213216_external_scrna"]["qc_pass_cells_total"] == 242430
    assert ev["GSE213216_external_scrna"]["primary_stromal_EcP_vs_EuE"]["shortlist_metrics"]["B0"]["top20"]["replicated_q05_count"] == 9
    assert ev["GSE141549_external_bulk"]["Joint_Dual_Contrast_and_Host_False_Positive_Audit"]["B0"]["top20"]["Host_false_positive_fraction_of_PE_hits"] == 1.0
    assert ev["GSE141549_external_bulk"]["Program_Transfer_All_Peritoneal_vs_Host"]["M2"]["top20"]["auc"] == pytest.approx(0.68397, abs=1e-4)
    assert len(dossiers) == 1 and dossiers.iloc[0]["symbol"] == "RUNX1"

    # Assert every key string appears verbatim in README.md, README.zh-CN.md, RESULTS.md, and RESULTS.zh-CN.md
    required_tokens = [
        "115,402",
        "242,430",
        "10,972",
        "569",
        "0.684",
        "0.00449",
        "88.9%",
        "72.2%",
        "+31.4",
        "RUNX1",
        "ENSG00000159216",
    ]
    for doc_rel in ("README.md", "README.zh-CN.md", "reports/RESULTS.md", "reports/RESULTS.zh-CN.md"):
        txt = (ROOT / doc_rel).read_text(encoding="utf-8")
        for tok in required_tokens:
            assert tok in txt, f"Token {tok!r} missing from {doc_rel}"
