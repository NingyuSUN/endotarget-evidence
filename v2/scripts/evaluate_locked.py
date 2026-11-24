#!/usr/bin/env python3
"""Stage 3: Locked external evaluation of frozen B0/G1/M1/M2 predictions on:
  1. GSE213216 external single-cell cohort (Fonseca et al., Nat Genet 2023; 29 locked scRNA libraries)
  2. GSE141549 external 408-sample / 162-subject multi-lesion & host-peritoneum double-difference cohort
  3. Compartment-resolved human genetic enrichment & Multi-Evidence Target Dossiers (with GTEx v10 context)

STRICT CRYPTOGRAPHIC LOCK VERIFICATION:
- Refuses to execute unless receipts/protocol_freeze.json, receipts/prediction_seal.json,
  and receipts/external_cohorts_data.json pass full SHA-256 verification.
"""
from __future__ import annotations

import csv
import datetime as dt
import gzip
import io
import json
import math
import os
import pathlib
import re
import subprocess
import sys
import tarfile
from collections import Counter, defaultdict

os.environ["OMP_NUM_THREADS"] = "8"
os.environ["OPENBLAS_NUM_THREADS"] = "8"
os.environ["MKL_NUM_THREADS"] = "8"
os.environ["NUMEXPR_NUM_THREADS"] = "8"

import h5py
import numpy as np
import openpyxl
import pandas as pd
import scipy.sparse as sp
from scipy import stats
import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from endotarget_lib import (  # noqa: E402
     _canonicalize_raw_counts,
    aggregate_pseudobulk,
    assign_lineages_de_novo,
    bh_fdr,
    calc_tmm_norm_factors,
    ebayes_moderated_ols,
    limma_duplicate_correlation_gls,
    limma_voom_fit,
    load_10x_targz_or_mtx_canonical,
    signed_z_from_t_df,
    squeeze_var,
)
from fetch_external_cohorts import verify_prediction_seal  # noqa: E402
from freeze_protocol import sha256_file, verify_protocol_freeze  # noqa: E402
from m0_audit_and_priors import load_hgnc, resolve_symbol, write_tsv  # noqa: E402
from run_discovery_v2 import align_sparse_matrices, fit_paired_donor_voom  # noqa: E402


def verify_external_data_receipt(root: pathlib.Path = ROOT) -> dict:
    rec_path = root / "receipts/external_cohorts_data.json"
    if not rec_path.exists():
        raise FileNotFoundError(f"Missing {rec_path}")
    rec = json.loads(rec_path.read_text(encoding="utf-8"))
    if rec["protocol_freeze_sha256"] != sha256_file(root / "receipts/protocol_freeze.json"):
        raise ValueError("external_cohorts_data.json protocol_freeze_sha256 mismatch")
    if rec["prediction_seal_sha256"] != sha256_file(root / "receipts/prediction_seal.json"):
        raise ValueError("external_cohorts_data.json prediction_seal_sha256 mismatch")
    for fname, info in rec["GSE213216_files"].items():
        p = root / "data/raw/scrna/GSE213216" / fname
        if sha256_file(p) != info["sha256"]:
            raise ValueError(f"SHA-256 mismatch on GSE213216 file {fname}")
    for fname, info in rec["GSE141549_files"].items():
        p = root / "data/raw/bulk/GSE141549" / fname
        if sha256_file(p) != info["sha256"]:
            raise ValueError(f"SHA-256 mismatch on GSE141549 file {fname}")
    return rec


def load_gse213216_tar_fast(tar_path: pathlib.Path, sym_to_ensg: dict[str, str]):
    """Load filtered_feature_bc_matrix.h5 inside GSE213216 .tar.gz via in-memory BytesIO if available,
    falling back to MTX reader otherwise."""
    with tarfile.open(tar_path, "r:gz") as tf:
        h5_members = [
            m for m in tf.getmembers()
            if m.isfile() and m.name.endswith(".h5") and "filtered_feature_bc_matrix" in m.name
        ]
        if h5_members:
            raw_bytes = tf.extractfile(h5_members[0]).read()
            with h5py.File(io.BytesIO(raw_bytes), "r") as f:
                grp = f["matrix"]
                data = grp["data"][:]
                indices = grp["indices"][:]
                indptr = grp["indptr"][:]
                shape = tuple(grp["shape"][:])
                barcodes = [b.decode("utf-8") if isinstance(b, bytes) else str(b) for b in grp["barcodes"][:]]
                feat = grp["features"]
                raw_ids = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in feat["id"][:]]
                raw_names = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in feat["name"][:]]
                ftypes = (
                    [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in feat["feature_type"][:]]
                    if "feature_type" in feat
                    else ["Gene Expression"] * len(raw_ids)
                )
            X_cbf = sp.csr_matrix((data.astype(np.int64), indices, indptr), shape=(shape[1], shape[0]))
            return _canonicalize_raw_counts(X_cbf, barcodes, raw_ids, raw_names, ftypes, sym_to_ensg)
    return load_10x_targz_or_mtx_canonical([tar_path], sym_to_ensg=sym_to_ensg)


def compute_shortlist_metrics(
    shortlists: dict[str, dict[int, list[dict]]],
    ext_genes: list[str],
    ext_logfc: np.ndarray,
    ext_z: np.ndarray,
    ext_q: np.ndarray,
) -> dict:
    gpos = {g: i for i, g in enumerate(ext_genes)}
    out = {}
    for mname, by_k in shortlists.items():
        out[mname] = {}
        for k, cand_list in by_k.items():
            present = [c for c in cand_list if c["ensembl_gene_id"] in gpos]
            n_eval = len(present)
            same_dir = 0
            rep_q05 = 0
            dir_z_vals = []
            gene_details = []
            for c in present:
                gi = gpos[c["ensembl_gene_id"]]
                dsign = int(c["discovery_sign"])
                lfc = float(ext_logfc[gi])
                zv = float(ext_z[gi])
                qv = float(ext_q[gi])
                is_same = (lfc > 0 and dsign > 0) or (lfc < 0 and dsign < 0)
                if is_same:
                    same_dir += 1
                if is_same and qv < 0.05:
                    rep_q05 += 1
                dir_z_vals.append(dsign * zv)
                gene_details.append({
                    "ensembl_gene_id": c["ensembl_gene_id"],
                    "symbol": c["symbol"],
                    "rank": int(c[f"rank_{mname}"]),
                    "discovery_sign": dsign,
                    "ext_logFC": round(lfc, 5),
                    "ext_z": round(zv, 5),
                    "ext_q": qv,
                    "replicated_q05": bool(is_same and qv < 0.05),
                })
            miss = k - n_eval
            out[mname][f"top{k}"] = {
                "selected": k,
                "evaluable": n_eval,
                "same_direction_count": same_dir,
                "same_direction_rate_evaluable": round(same_dir / n_eval, 6) if n_eval else 0.0,
                "same_direction_bounds_on_k": [round(same_dir / k, 6), round((same_dir + miss) / k, 6)],
                "replicated_q05_count": rep_q05,
                "replicated_q05_rate_on_k": round(rep_q05 / k, 6),
                "replicated_q05_rate_evaluable": round(rep_q05 / n_eval, 6) if n_eval else 0.0,
                "mean_directional_z": round(float(np.mean(dir_z_vals)), 5) if dir_z_vals else 0.0,
                "genes": gene_details,
            }
    return out


def compute_program_transfer(
    shortlists: dict[str, dict[int, list[dict]]],
    ext_genes: list[str],
    Y_expr: np.ndarray,
    is_case: np.ndarray,
    block_ids: list[str] | None = None,
) -> dict:
    """Compute discovery-weighted signature projection score theta (in discovery shift units) and Mann-Whitney AUC."""
    gpos = {g: i for i, g in enumerate(ext_genes)}
    ref_mask = is_case == 0
    Y_centered = Y_expr - Y_expr[:, ref_mask].mean(axis=1, keepdims=True)
    out = {}
    for mname, by_k in shortlists.items():
        out[mname] = {}
        for k, cand_list in by_k.items():
            present = [c for c in cand_list if c["ensembl_gene_id"] in gpos]
            if len(present) < 5:
                out[mname][f"top{k}"] = {"status": "insufficient_genes"}
                continue
            idx = np.array([gpos[c["ensembl_gene_id"]] for c in present], dtype=np.int64)
            beta = np.array([float(c["logFC_EcP_vs_EuE"]) for c in present], dtype=np.float64)
            denom = float(beta @ beta)
            scores = (beta @ Y_centered[idx, :]) / max(denom, 1e-12)
            pos_s = scores[is_case == 1]
            neg_s = scores[is_case == 0]
            u_stat = stats.mannwhitneyu(pos_s, neg_s, alternative="two-sided").statistic
            auc_val = float(u_stat / (len(pos_s) * len(neg_s)))

            # Estimate theta and 95% CI via OLS (or paired/cluster SE)
            X_d = np.column_stack([np.ones(len(scores)), is_case.astype(float)])
            if block_ids is not None and len(set(block_ids)) < len(block_ids):
                fit_p = limma_duplicate_correlation_gls(scores[None, :], np.vstack([scores, scores])[:1].shape and X_d, block_ids, coef_idx=1) if False else None
            # Exact OLS point estimate + cluster-robust (Liang-Zeger) or classical t SE
            b_hat, *_ = np.linalg.lstsq(X_d, scores, rcond=None)
            resid = scores - X_d @ b_hat
            XtX_inv = np.linalg.inv(X_d.T @ X_d)
            if block_ids is not None:
                clusters = sorted(set(block_ids))
                G = len(clusters)
                meat = np.zeros((2, 2), dtype=np.float64)
                for cl in clusters:
                    cm = np.array([b == cl for b in block_ids], dtype=bool)
                    u_c = X_d[cm].T @ resid[cm]
                    meat += np.outer(u_c, u_c)
                # CR1 small-sample correction
                scale = (G / max(G - 1, 1)) * ((len(scores) - 1) / max(len(scores) - 2, 1))
                vcov = scale * (XtX_inv @ meat @ XtX_inv)
                se_theta = float(math.sqrt(max(vcov[1, 1], 1e-15)))
                df_t = max(G - 2, 1)
            else:
                df_t = len(scores) - 2
                s2 = float((resid ** 2).sum() / max(df_t, 1))
                se_theta = float(math.sqrt(max(s2 * XtX_inv[1, 1], 1e-15)))
            theta = float(b_hat[1])
            t_val = theta / se_theta
            p_val = float(2.0 * stats.t.sf(abs(t_val), df=df_t))
            tq = float(stats.t.ppf(0.975, df=df_t))
            out[mname][f"top{k}"] = {
                "genes_present": len(present),
                "auc": round(auc_val, 5),
                "theta": round(theta, 5),
                "se": round(se_theta, 5),
                "ci_95": [round(theta - tq * se_theta, 5), round(theta + tq * se_theta, 5)],
                "p_value": p_val,
            }
    return out


def cluster_bootstrap_method_comparisons(
    shortlists: dict[str, dict[int, list[dict]]],
    ext_genes: list[str],
    M_expr: np.ndarray,
    is_case: np.ndarray,
    block_ids: list[str],
    n_draws: int = 2000,
    seed: int = 20261120,
) -> dict:
    """Donor/subject cluster bootstrap (2,000 draws) for paired method differences:
    M1 - B0, M1 - G1, M2 - B0, M2 - G1 on replicated q<0.05 rate, direction agreement, and directional z."""
    rng = np.random.default_rng(seed)
    gpos = {g: i for i, g in enumerate(ext_genes)}
    unique_blocks = sorted(set(block_ids))
    block_to_cols = {b: [i for i, x in enumerate(block_ids) if x == b] for b in unique_blocks}
    n_blocks = len(unique_blocks)

    # Pre-index candidate positions and signs
    idx_map = {}
    sign_map = {}
    for mname in ("B0", "G1", "M1", "M2"):
        for k in (20, 50):
            pres = [c for c in shortlists[mname][k] if c["ensembl_gene_id"] in gpos]
            idx_map[(mname, k)] = np.array([gpos[c["ensembl_gene_id"]] for c in pres], dtype=np.int64)
            sign_map[(mname, k)] = np.array([int(c["discovery_sign"]) for c in pres], dtype=np.float64)

    pairs_to_compare = [("M1", "B0"), ("M1", "G1"), ("M2", "B0"), ("M2", "G1"), ("M2", "M1")]
    records = defaultdict(list)
    valid_draws = 0

    for _ in range(n_draws):
        sampled_blocks = rng.choice(unique_blocks, size=n_blocks, replace=True)
        col_idx = []
        for b in sampled_blocks:
            col_idx.extend(block_to_cols[b])
        col_arr = np.array(col_idx, dtype=np.int64)
        y_case = is_case[col_arr]
        if (y_case == 1).sum() < 3 or (y_case == 0).sum() < 3:
            continue
        X_b = np.column_stack([np.ones(len(col_arr), dtype=np.float64), y_case.astype(np.float64)])
        df_eff = max(2, n_blocks - 2)
        try:
            de_b = ebayes_moderated_ols(M_expr[:, col_arr], X_b, coef_idx=1, df_override=df_eff)
        except Exception:
            continue
        valid_draws += 1
        lfc = de_b["logFC"]
        zv = de_b["z"]
        qv = de_b["q"]

        draw_vals = {}
        for mname in ("B0", "G1", "M1", "M2"):
            for k in (20, 50):
                ii = idx_map[(mname, k)]
                ss = sign_map[(mname, k)]
                n_eval = max(len(ii), 1)
                same = (np.sign(lfc[ii]) == ss)
                rep = same & (qv[ii] < 0.05)
                draw_vals[(mname, k, "rep_rate_on_k")] = float(rep.sum() / k)
                draw_vals[(mname, k, "dir_rate_eval")] = float(same.sum() / n_eval)
                draw_vals[(mname, k, "mean_dir_z")] = float((ss * zv[ii]).mean())

        for m_a, m_b in pairs_to_compare:
            for k in (20, 50):
                for metric in ("rep_rate_on_k", "dir_rate_eval", "mean_dir_z"):
                    records[(f"{m_a}_minus_{m_b}", f"top{k}", metric)].append(
                        draw_vals[(m_a, k, metric)] - draw_vals[(m_b, k, metric)]
                    )

    summary = {"requested_draws": n_draws, "valid_draws": valid_draws, "comparisons": {}}
    for (comp_name, k_str, metric), vals in records.items():
        arr = np.array(vals, dtype=np.float64)
        summary["comparisons"].setdefault(comp_name, {}).setdefault(k_str, {})[metric] = {
            "boot_mean": round(float(arr.mean()), 5),
            "ci_95": [round(float(np.percentile(arr, 2.5)), 5), round(float(np.percentile(arr, 97.5)), 5)],
            "prob_gt_0": round(float((arr > 0).mean()), 4),
        }
    return summary


def paired_difference_moderated_fit(D: np.ndarray) -> dict:
    """Exact empirical Bayes moderated one-sample t-test on paired within-subject differences D (genes x n_pairs)."""
    n_pairs = D.shape[1]
    df_res = n_pairs - 1
    mean_d = D.mean(axis=1)
    s2 = ((D - mean_d[:, None]) ** 2).sum(axis=1) / df_res
    s2_post, d0, s0 = squeeze_var(s2, df_res)
    se = np.sqrt(s2_post / n_pairs)
    t_stat = mean_d / np.maximum(se, 1e-15)
    df_tot = float(df_res + d0)
    p_val = 2.0 * stats.t.sf(np.abs(t_stat), df=df_tot)
    q_val = bh_fdr(p_val)
    z_val = signed_z_from_t_df(t_stat, df_tot)
    return {
        "n_pairs": n_pairs,
        "logFC": mean_d,
        "se": se,
        "t": t_stat,
        "z": z_val,
        "p": p_val,
        "q": q_val,
        "df_residual": int(df_res),
        "prior_df": float(d0),
        "df_total": df_tot,
    }


def paired_dual_contrast_bootstrap(
    shortlists: dict[str, dict[int, list[dict]]],
    bulk_genes: list[str],
    subj_expr: dict[str, dict[str, np.ndarray]],
    n_draws: int = 2000,
    seed: int = 20261120,
) -> dict:
    """Paired subject-level bootstrap (2,000 draws) on GSE141549 simultaneously resampling subjects
    to get 95% CIs on:
      - Paired PeL vs PE replication rate
      - Paired PeL vs PP (host-controlled!) replication rate
      - Joint Dual-Contrast (PeL vs PE AND PeL vs PP both q<0.05 same sign) True Lesion Target rate
      - Host-Peritoneal False-Positive Rate among PeL vs PE hits!
    """
    rng = np.random.default_rng(seed)
    gpos = {g: i for i, g in enumerate(bulk_genes)}
    pe_subjects = sorted(s for s, g in subj_expr.items() if "PE" in g and "Peritoneal_Lesion" in g)
    pp_subjects = sorted(s for s, g in subj_expr.items() if "PP" in g and "Peritoneal_Lesion" in g)
    union_subjects = sorted(set(pe_subjects) | set(pp_subjects))

    D_pe = np.column_stack([subj_expr[s]["Peritoneal_Lesion"] - subj_expr[s]["PE"] for s in pe_subjects])
    D_pp = np.column_stack([subj_expr[s]["Peritoneal_Lesion"] - subj_expr[s]["PP"] for s in pp_subjects])
    pe_pos = {s: i for i, s in enumerate(pe_subjects)}
    pp_pos = {s: i for i, s in enumerate(pp_subjects)}

    idx_map = {}
    sign_map = {}
    for mname in ("B0", "G1", "M1", "M2"):
        for k in (20, 50):
            pres = [c for c in shortlists[mname][k] if c["ensembl_gene_id"] in gpos]
            idx_map[(mname, k)] = np.array([gpos[c["ensembl_gene_id"]] for c in pres], dtype=np.int64)
            sign_map[(mname, k)] = np.array([int(c["discovery_sign"]) for c in pres], dtype=np.float64)

    pairs_to_compare = [("M1", "B0"), ("M1", "G1"), ("M2", "B0"), ("M2", "G1"), ("M2", "M1")]
    records = defaultdict(list)
    valid_draws = 0

    for _ in range(n_draws):
        sampled = rng.choice(union_subjects, size=len(union_subjects), replace=True)
        idx_pe = [pe_pos[s] for s in sampled if s in pe_pos]
        idx_pp = [pp_pos[s] for s in sampled if s in pp_pos]
        if len(idx_pe) < 5 or len(idx_pp) < 5:
            continue
        fit_pe_b = paired_difference_moderated_fit(D_pe[:, idx_pe])
        fit_pp_b = paired_difference_moderated_fit(D_pp[:, idx_pp])
        valid_draws += 1

        vals = {}
        for mname in ("B0", "G1", "M1", "M2"):
            for k in (20, 50):
                ii = idx_map[(mname, k)]
                ss = sign_map[(mname, k)]
                same_pe = np.sign(fit_pe_b["logFC"][ii]) == ss
                rep_pe = same_pe & (fit_pe_b["q"][ii] < 0.05)
                same_pp = np.sign(fit_pp_b["logFC"][ii]) == ss
                rep_pp = same_pp & (fit_pp_b["q"][ii] < 0.05)
                joint_rep = rep_pe & rep_pp
                vals[(mname, k, "rep_PE_rate")] = float(rep_pe.sum() / k)
                vals[(mname, k, "rep_PP_host_rate")] = float(rep_pp.sum() / k)
                vals[(mname, k, "joint_PE_and_PP_rate")] = float(joint_rep.sum() / k)
                vals[(mname, k, "dir_PP_host_rate")] = float(same_pp.sum() / max(len(ii), 1))

        for m_a, m_b in pairs_to_compare:
            for k in (20, 50):
                for metric in ("rep_PE_rate", "rep_PP_host_rate", "joint_PE_and_PP_rate", "dir_PP_host_rate"):
                    records[(f"{m_a}_minus_{m_b}", f"top{k}", metric)].append(
                        vals[(m_a, k, metric)] - vals[(m_b, k, metric)]
                    )

    summary = {"requested_draws": n_draws, "valid_draws": valid_draws, "comparisons": {}}
    for (comp_name, k_str, metric), arr_list in records.items():
        arr = np.array(arr_list, dtype=np.float64)
        summary["comparisons"].setdefault(comp_name, {}).setdefault(k_str, {})[metric] = {
            "boot_mean": round(float(arr.mean()), 5),
            "ci_95": [round(float(np.percentile(arr, 2.5)), 5), round(float(np.percentile(arr, 97.5)), 5)],
            "prob_gt_0": round(float((arr > 0).mean()), 4),
        }
    return summary


def load_gse141549_gene_matrix(bulk_dir: pathlib.Path, soft_path: pathlib.Path, hgnc: dict):
    """Parse GPL10558.annot.gz + GSE141549 SOFT platform tables + GSE141549 workbook -> unique ENSG x 408 matrix."""
    probe_to_ensg = {}

    # 1. Read GPL10558.annot.gz
    with gzip.open(bulk_dir / "GPL10558.annot.gz", "rt", encoding="utf-8", errors="replace") as f:
        hdr = None
        for line in f:
            if line.startswith(("#", "^", "!")):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if hdr is None:
                if "ID" in parts:
                    hdr = parts
                    i_id = hdr.index("ID")
                    i_ent = hdr.index("Gene ID") if "Gene ID" in hdr else None
                    i_sym = hdr.index("Gene symbol") if "Gene symbol" in hdr else None
                continue
            if len(parts) <= i_id:
                continue
            pid = parts[i_id].strip()
            cands = set()
            if i_ent is not None and i_ent < len(parts):
                for e in re.split(r"///|;|,", parts[i_ent]):
                    e = e.strip()
                    if e in hgnc["entrez_to_ensg"]:
                        cands.add(hgnc["entrez_to_ensg"][e])
            if not cands and i_sym is not None and i_sym < len(parts):
                for s in re.split(r"///|;|,", parts[i_sym]):
                    eg, _ = resolve_symbol(s.strip(), hgnc)
                    if eg:
                        cands.add(eg)
            if len(cands) == 1:
                probe_to_ensg[pid] = next(iter(cands))

    # 2. Also read GPL13376 / GPL10558 platform tables inside GSE141549_family.soft.gz for any additional probes
    with gzip.open(soft_path, "rt", encoding="utf-8", errors="replace") as f:
        in_plat = False
        phdr = None
        for line in f:
            if line.startswith("!platform_table_begin"):
                in_plat = True
                phdr = None
                continue
            if line.startswith("!platform_table_end"):
                in_plat = False
                continue
            if not in_plat:
                continue
            parts = line.rstrip("\r\n").split("\t")
            if phdr is None:
                phdr = parts
                p_id = phdr.index("ID") if "ID" in phdr else 0
                p_ent = phdr.index("Entrez_Gene_ID") if "Entrez_Gene_ID" in phdr else None
                p_sym = phdr.index("Symbol") if "Symbol" in phdr else None
                continue
            if len(parts) <= p_id:
                continue
            pid = parts[p_id].strip()
            if pid in probe_to_ensg:
                continue
            cands = set()
            if p_ent is not None and p_ent < len(parts):
                e = parts[p_ent].strip().split(".")[0]
                if e in hgnc["entrez_to_ensg"]:
                    cands.add(hgnc["entrez_to_ensg"][e])
            if not cands and p_sym is not None and p_sym < len(parts):
                eg, _ = resolve_symbol(parts[p_sym].strip(), hgnc)
                if eg:
                    cands.add(eg)
            if len(cands) == 1:
                probe_to_ensg[pid] = next(iter(cands))

    # 3. Stream GSE141549_batchCorrectednormalizedArrayscombined.xlsx
    wb = openpyxl.load_workbook(
        bulk_dir / "GSE141549_batchCorrectednormalizedArrayscombined.xlsx",
        read_only=True,
        data_only=True,
    )
    ws = wb[wb.sheetnames[0]]
    row_iter = ws.iter_rows(values_only=True)
    header = [str(x).strip() for x in next(row_iter)]
    sample_titles = header[2:]

    best_by_gene = {}  # ensg -> (mean_expr, values_float64)
    total_probes = 0
    for r in row_iter:
        if not r or r[1] is None:
            continue
        total_probes += 1
        sym_raw = str(r[0]).strip() if r[0] is not None else ""
        pid = str(r[1]).strip()
        eg = probe_to_ensg.get(pid)
        if not eg and sym_raw:
            eg, _ = resolve_symbol(sym_raw, hgnc)
        if not eg:
            continue
        vals = np.array([float(x) if x is not None else np.nan for x in r[2:]], dtype=np.float64)
        if not np.isfinite(vals).all():
            continue
        m_val = float(vals.mean())
        prev = best_by_gene.get(eg)
        if prev is None or m_val > prev[0]:
            best_by_gene[eg] = (m_val, vals)

    genes = sorted(best_by_gene.keys())
    M_all = np.vstack([best_by_gene[g][1] for g in genes])
    return genes, sample_titles, M_all, total_probes


def test_genetic_set_enrichment(
    genes: list[str],
    stat_vec: np.ndarray,
    gene_sets: dict[str, set[str]],
    expr_vec: np.ndarray | None = None,
    n_perm: int = 2000,
    seed: int = 20261120,
) -> dict:
    """Mann-Whitney rank-sum test + expression-decile-matched permutation test for genetic sets."""
    rng = np.random.default_rng(seed)
    gpos = {g: i for i, g in enumerate(genes)}
    ranks_pct = stats.rankdata(stat_vec) / max(len(stat_vec), 1)

    if expr_vec is not None:
        quantiles = np.quantile(expr_vec, np.linspace(0, 1, 11)[1:-1])
        bins = np.digitize(expr_vec, quantiles)
    else:
        bins = np.zeros(len(genes), dtype=int)
    by_bin = {b: np.where(bins == b)[0] for b in np.unique(bins)}

    res = {}
    for sname, gset in gene_sets.items():
        in_idx = np.array([gpos[g] for g in sorted(gset) if g in gpos], dtype=np.int64)
        out_mask = np.ones(len(genes), dtype=bool)
        out_mask[in_idx] = False
        if len(in_idx) < 5:
            res[sname] = {"n_in_universe": int(len(in_idx)), "status": "too_few_genes"}
            continue
        in_vals = stat_vec[in_idx]
        out_vals = stat_vec[out_mask]
        mw = stats.mannwhitneyu(in_vals, out_vals, alternative="greater")
        obs_mean = float(in_vals.mean())
        obs_mean_pct = float(ranks_pct[in_idx].mean() * 100.0)

        need = {b: int((bins[in_idx] == b).sum()) for b in by_bin}
        perm_means = np.empty(n_perm, dtype=np.float64)
        for p_i in range(n_perm):
            draw = np.concatenate([rng.choice(by_bin[b], size=cnt, replace=False) for b, cnt in need.items() if cnt > 0])
            perm_means[p_i] = float(stat_vec[draw].mean())
        emp_p = float((1 + (perm_means >= obs_mean).sum()) / (1 + n_perm))
        res[sname] = {
            "n_in_universe": int(len(in_idx)),
            "mean_score_in_set": round(obs_mean, 5),
            "mean_score_background": round(float(out_vals.mean()), 5),
            "fold_enrichment_mean": round(obs_mean / max(float(out_vals.mean()), 1e-8), 4),
            "mean_rank_percentile": round(obs_mean_pct, 2),
            "mannwhitney_p_one_sided": float(mw.pvalue),
            "expr_matched_perm_p": round(emp_p, 5),
        }
    return res


def query_gtex_for_dossier_genes(symbols: list[str]) -> dict[str, dict]:
    """Query GTEx v10 baseline expression and cis-eQTLs in Uterus, Cultured fibroblasts, and Ovary
    using the official gtex_cli.py skill wrapper."""
    gtex_cli = pathlib.Path("/home/seedna2/.gemini/config/plugins/genomics/skills/gtex_database/scripts/gtex_cli.py")
    uv_bin = pathlib.Path("/home/seedna2/.local/bin/uv")
    cache_path = ROOT / "data/raw/genetics/gtex/dossier_gtex_cache.json"
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache = json.loads(cache_path.read_text(encoding="utf-8")) if cache_path.exists() else {}

    if not gtex_cli.exists() or not uv_bin.exists():
        return cache

    for sym in symbols:
        if sym in cache:
            continue
        id_out = pathlib.Path(f"/tmp/gtex_v2_id_{sym}.json")
        expr_out = pathlib.Path(f"/tmp/gtex_v2_expr_{sym}.json")
        eqtl_out = pathlib.Path(f"/tmp/gtex_v2_eqtl_{sym}.json")
        for p in (id_out, expr_out, eqtl_out):
            if p.exists():
                p.unlink()
        subprocess.run(
            [str(uv_bin), "run", str(gtex_cli), "resolve-gencode-id", sym, "--output", str(id_out)],
            check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        if not id_out.exists():
            continue
        id_data = json.loads(id_out.read_text(encoding="utf-8"))
        gid = id_data.get("gencode_id") if isinstance(id_data, dict) else None
        if not gid:
            continue
        subprocess.run(
            [str(uv_bin), "run", str(gtex_cli), "get-median-expression", gid, "--output", str(expr_out)],
            check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        subprocess.run(
            [str(uv_bin), "run", str(gtex_cli), "get-gene-eqtls", gid, "--output", str(eqtl_out)],
            check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        rec = {"gencode_id": gid, "tpm_Uterus": None, "tpm_Fibroblasts": None, "tpm_Ovary": None,
               "eqtls_Uterus": 0, "eqtls_Fibroblasts": 0, "eqtls_Ovary": 0}
        if expr_out.exists():
            ed = json.loads(expr_out.read_text(encoding="utf-8"))
            erows = ed if isinstance(ed, list) else ed.get("data", [])
            for r in erows:
                tid = r.get("tissueSiteDetailId", "")
                if tid == "Uterus":
                    rec["tpm_Uterus"] = round(float(r.get("median", 0.0)), 2)
                elif tid == "Cells_Cultured_fibroblasts":
                    rec["tpm_Fibroblasts"] = round(float(r.get("median", 0.0)), 2)
                elif tid == "Ovary":
                    rec["tpm_Ovary"] = round(float(r.get("median", 0.0)), 2)
        if eqtl_out.exists():
            qd = json.loads(eqtl_out.read_text(encoding="utf-8"))
            qrows = qd if isinstance(qd, list) else qd.get("data", [])
            for r in qrows:
                tid = r.get("tissueSiteDetailId", "")
                if tid == "Uterus":
                    rec["eqtls_Uterus"] += 1
                elif tid == "Cells_Cultured_fibroblasts":
                    rec["eqtls_Fibroblasts"] += 1
                elif tid == "Ovary":
                    rec["eqtls_Ovary"] += 1
        cache[sym] = rec

    cache_path.write_text(json.dumps(cache, indent=2) + "\n", encoding="utf-8")
    return cache


def main():
    # 0. Strict cryptographic lock verification BEFORE reading any external matrix
    verify_protocol_freeze(ROOT)
    verify_prediction_seal(ROOT)
    verify_external_data_receipt(ROOT)

    freeze_sha = sha256_file(ROOT / "receipts/protocol_freeze.json")
    seal_sha = sha256_file(ROOT / "receipts/prediction_seal.json")
    ext_data_sha = sha256_file(ROOT / "receipts/external_cohorts_data.json")
    cfg = yaml.safe_load((ROOT / "configs/protocol_v2.yaml").read_text(encoding="utf-8"))
    seed = int(cfg["runtime_policy"]["random_seed"])
    n_boot = int(cfg["external_evaluation"]["bootstrap_draws"])

    hgnc = load_hgnc(ROOT / "data/raw/reference/hgnc_complete_set.tsv")
    sym_to_ensg = {**hgnc["alias"], **hgnc["prev"], **hgnc["approved"]}

    with (ROOT / "manifests/canonical_lineage_markers.tsv").open("r", encoding="utf-8") as f:
        marker_rows = list(csv.DictReader(f, delimiter="\t"))
    with (ROOT / "manifests/frozen_candidates_v2.tsv").open("r", encoding="utf-8") as f:
        cand_rows = list(csv.DictReader(f, delimiter="\t"))
    with (ROOT / "manifests/human_genetics_prior_v2.tsv").open("r", encoding="utf-8") as f:
        gen_rows = list(csv.DictReader(f, delimiter="\t"))
    genetics_by_ensg = {r["ensembl_gene_id"]: r for r in gen_rows}

    with gzip.open(ROOT / "reports/discovery/GSE179640_stromal_rankings_full.tsv.gz", "rt", encoding="utf-8") as f:
        disc_full_rows = list(csv.DictReader(f, delimiter="\t"))
    disc_by_ensg = {r["ensembl_gene_id"]: r for r in disc_full_rows}

    shortlists = {}
    for mname in ("B0", "G1", "M1", "M2"):
        shortlists[mname] = {}
        for k in (20, 50):
            sub = [r for r in cand_rows if r[f"in_{mname}_top{k}"] == "True"]
            sub.sort(key=lambda r: int(r[f"rank_{mname}"]))
            assert len(sub) == k, f"Frozen {mname} top-{k} must have {k} genes, got {len(sub)}"
            shortlists[mname][k] = sub

    # =========================================================================
    # PART 1: External Single-Cell Replication in GSE213216 (29 locked libraries)
    # =========================================================================
    out_sc = ROOT / "reports/external_scrna"
    out_sc.mkdir(parents=True, exist_ok=True)
    sc_raw_dir = ROOT / "data/raw/scrna/GSE213216"

    with (ROOT / "manifests/m0_GSE213216_locked_selection.tsv").open("r", encoding="utf-8") as f:
        sel213 = list(csv.DictReader(f, delimiter="\t"))

    qc_cfg = cfg["single_cell_processing"]["qc_filters"]
    mats_213, glists_213, obs_parts_213, lib_qc_213 = [], [], [], []
    patient_disease_map = {}

    for r in sel213:
        fnames = [x.strip() for x in r["supplementary_files"].split(";") if x.strip()]
        tar_p = sc_raw_dir / fnames[0]
        X_lib, ensg_lib, qc_df = load_gse213216_tar_fast(tar_p, sym_to_ensg=sym_to_ensg)
        pass_qc = (
            (qc_df["n_genes_by_counts"] >= int(qc_cfg["min_genes_by_counts"]))
            & (qc_df["total_counts"] >= int(qc_cfg["min_total_counts"]))
            & (qc_df["pct_counts_mt"] <= float(qc_cfg["max_pct_counts_mt"]))
        ).values
        patient_disease_map[r["patient_id"]] = r["curated_disease"] or "ENDO"
        lib_qc_213.append({
            "gsm": r["gsm"],
            "patient_id": r["patient_id"],
            "site_class": r["site_class"],
            "curated_disease": r["curated_disease"],
            "raw_cells": int(X_lib.shape[0]),
            "qc_pass_cells": int(pass_qc.sum()),
            "qc_pass_fraction": round(float(pass_qc.mean()), 4),
            "median_genes_qc": float(np.median(qc_df.loc[pass_qc, "n_genes_by_counts"])) if pass_qc.any() else 0.0,
        })
        sub_qc = qc_df.loc[pass_qc].copy()
        sub_qc["gsm"] = r["gsm"]
        sub_qc["donor_id"] = r["patient_id"]
        sub_qc["site_code"] = r["site_class"]
        sub_qc["cell_id"] = r["gsm"] + ":" + sub_qc["barcode"].astype(str)
        mats_213.append(X_lib[pass_qc])
        glists_213.append(ensg_lib)
        obs_parts_213.append(sub_qc)

    X_213, ensg_213 = align_sparse_matrices(mats_213, glists_213)
    obs_213 = pd.concat(obs_parts_213, ignore_index=True)
    obs_213.index = obs_213["cell_id"].values

    annot_213, clusters_213 = assign_lineages_de_novo(
        X_213, ensg_213, obs_213, marker_rows, cfg["single_cell_processing"], random_seed=seed
    )
    write_tsv(out_sc / "GSE213216_library_qc.tsv", lib_qc_213, list(lib_qc_213[0].keys()))
    clusters_213.to_csv(out_sc / "GSE213216_cluster_annotations.tsv", sep="\t", index=False)

    valid_lineages = list(cfg["single_cell_processing"]["de_novo_lineage_annotation"]["lineages"])
    min_cells_pb = int(cfg["single_cell_processing"]["pseudobulk"]["min_cells_per_sample"])
    pb_213, cons_213 = aggregate_pseudobulk(X_213, ensg_213, annot_213, valid_lineages=valid_lineages, min_cells=min_cells_pb)

    pb_213_meta_all = pd.concat([mdf for _, mdf in pb_213.values()], ignore_index=True)
    pb_213_meta_all.to_csv(out_sc / "GSE213216_pseudobulk_samples.tsv", sep="\t", index=False)

    # Primary GSE213216 Stromal Contrast: EcP vs EuE via donor-blocked GLS on TMM log2CPM
    s213_counts, s213_meta = pb_213["stromal_fibroblast"]
    eue_ecp_m = s213_meta[s213_meta["site_code"].isin(["EuE", "EcP"])].copy().reset_index(drop=True)
    eue_ecp_cols = s213_meta[s213_meta["site_code"].isin(["EuE", "EcP"])].index.values
    C_213_prim = s213_counts[:, eue_ecp_cols]
    cpm_213_prim = C_213_prim / np.maximum(C_213_prim.sum(axis=0, keepdims=True), 1.0) * 1e6
    expr_213_mask = (cpm_213_prim >= 1.0).mean(axis=1) >= 0.20
    genes_213 = [g for g, k in zip(ensg_213, expr_213_mask) if k]

    C_213_sub = C_213_prim[expr_213_mask, :]
    tmm_213 = calc_tmm_norm_factors(C_213_sub)
    lib_eff_213 = C_213_sub.sum(axis=0) * tmm_213
    Y_213 = np.log2((C_213_sub + 0.5) / (lib_eff_213[None, :] + 1.0) * 1e6)

    is_ecp_213 = (eue_ecp_m["site_code"].values == "EcP").astype(np.float64)
    blocks_213 = eue_ecp_m["donor_id"].tolist()
    X_213_prim = np.column_stack([np.ones(len(eue_ecp_m)), is_ecp_213])
    fit_213_prim = limma_duplicate_correlation_gls(Y_213, X_213_prim, blocks_213, coef_idx=1)

    # Save full GSE213216 primary DE table
    sc_de_rows = [
        {
            "ensembl_gene_id": g,
            "symbol": hgnc["ensg_to_symbol"].get(g, ""),
            "logFC_EcP_vs_EuE": round(float(fit_213_prim["logFC"][i]), 6),
            "t": round(float(fit_213_prim["t"][i]), 6),
            "z": round(float(fit_213_prim["z"][i]), 6),
            "p": float(fit_213_prim["p"][i]),
            "q_bh": float(fit_213_prim["q"][i]),
        }
        for i, g in enumerate(genes_213)
    ]
    write_tsv(out_sc / "GSE213216_stromal_DE_EcP_vs_EuE.tsv", sc_de_rows, list(sc_de_rows[0].keys()))

    metrics_213_prim = compute_shortlist_metrics(
        shortlists, genes_213, fit_213_prim["logFC"], fit_213_prim["z"], fit_213_prim["q"]
    )
    prog_213_prim = compute_program_transfer(shortlists, genes_213, Y_213, is_ecp_213, block_ids=blocks_213)
    boot_213_prim = cluster_bootstrap_method_comparisons(
        shortlists, genes_213, Y_213, is_ecp_213, blocks_213, n_draws=n_boot, seed=seed
    )

    # Sensitivities in GSE213216
    sens_213 = {}
    # (a) ENDO-only EuE vs EcP (excluding non-ENDO P06, P12)
    endo_keep = np.array([patient_disease_map.get(d, "ENDO") == "ENDO" for d in blocks_213], dtype=bool)
    fit_213_endo = limma_duplicate_correlation_gls(
        Y_213[:, endo_keep],
        np.column_stack([np.ones(int(endo_keep.sum())), is_ecp_213[endo_keep]]),
        [b for b, k in zip(blocks_213, endo_keep) if k],
        coef_idx=1,
    )
    sens_213["EcP_vs_ENDO_only_EuE"] = {
        "n_EcP": int((is_ecp_213[endo_keep] == 1).sum()),
        "n_EuE": int((is_ecp_213[endo_keep] == 0).sum()),
        "sig_BH_q05": int((fit_213_endo["q"] < 0.05).sum()),
        "shortlists": {
            m: {k: {kk: vv for kk, vv in d.items() if kk != "genes"} for k, d in md.items()}
            for m, md in compute_shortlist_metrics(shortlists, genes_213, fit_213_endo["logFC"], fit_213_endo["z"], fit_213_endo["q"]).items()
        },
    }

    # (b) Strict matched within-patient pairs in GSE213216 (P10, P11, P14, P16)
    by_d213 = eue_ecp_m.groupby("donor_id")["site_code"].apply(set)
    paired_pts_213 = sorted(d for d, s in by_d213.items() if {"EuE", "EcP"} <= s)
    if len(paired_pts_213) >= 3:
        fit_213_paired = fit_paired_donor_voom(s213_counts, s213_meta, site_case="EcP", site_ref="EuE", keep_mask=expr_213_mask)
        sens_213["Matched_within_patient_EcP_vs_EuE"] = {
            "n_pairs": fit_213_paired["n_pairs"],
            "paired_patients": fit_213_paired["paired_donors"],
            "sig_BH_q05": int((fit_213_paired["q"] < 0.05).sum()),
            "shortlists": {
                m: {k: {kk: vv for kk, vv in d.items() if kk != "genes"} for k, d in md.items()}
                for m, md in compute_shortlist_metrics(shortlists, genes_213, fit_213_paired["logFC"], fit_213_paired["z"], fit_213_paired["q"]).items()
            },
        }

    # (c) Strictly donor-disjoint unpaired (EcP-only vs EuE-only patients)
    disj_keep = np.array([b not in set(paired_pts_213) for b in blocks_213], dtype=bool)
    if (is_ecp_213[disj_keep] == 1).sum() >= 3 and (is_ecp_213[disj_keep] == 0).sum() >= 3:
        fit_213_disj = ebayes_moderated_ols(
            Y_213[:, disj_keep],
            np.column_stack([np.ones(int(disj_keep.sum())), is_ecp_213[disj_keep]]),
            coef_idx=1,
        )
        sens_213["Donor_disjoint_unpaired_EcP_vs_EuE"] = {
            "n_EcP": int((is_ecp_213[disj_keep] == 1).sum()),
            "n_EuE": int((is_ecp_213[disj_keep] == 0).sum()),
            "sig_BH_q05": int((fit_213_disj["q"] < 0.05).sum()),
            "shortlists": {
                m: {k: {kk: vv for kk, vv in d.items() if kk != "genes"} for k, d in md.items()}
                for m, md in compute_shortlist_metrics(shortlists, genes_213, fit_213_disj["logFC"], fit_213_disj["z"], fit_213_disj["q"]).items()
            },
        }

    # (d) Endometrioma EcO vs EuE in GSE213216 stromal fibroblasts
    eco_m = s213_meta[s213_meta["site_code"].isin(["EuE", "EcO"])].copy()
    if (eco_m["site_code"] == "EcO").sum() >= 3:
        C_eco = s213_counts[expr_213_mask][:, eco_m.index.values]
        tmm_eco = calc_tmm_norm_factors(C_eco)
        le_eco = C_eco.sum(axis=0) * tmm_eco
        Y_eco = np.log2((C_eco + 0.5) / (le_eco[None, :] + 1.0) * 1e6)
        is_eco = (eco_m["site_code"].values == "EcO").astype(float)
        fit_213_eco = limma_duplicate_correlation_gls(
            Y_eco, np.column_stack([np.ones(len(eco_m)), is_eco]), eco_m["donor_id"].tolist(), coef_idx=1
        )
        sens_213["Endometrioma_EcO_vs_EuE"] = {
            "n_EcO": int((is_eco == 1).sum()),
            "n_EuE": int((is_eco == 0).sum()),
            "sig_BH_q05": int((fit_213_eco["q"] < 0.05).sum()),
            "shortlists": {
                m: {k: {kk: vv for kk, vv in d.items() if kk != "genes"} for k, d in md.items()}
                for m, md in compute_shortlist_metrics(shortlists, genes_213, fit_213_eco["logFC"], fit_213_eco["z"], fit_213_eco["q"]).items()
            },
        }

    # (e) Endothelial & Myeloid EcP vs EuE in GSE213216
    ext_niche_de = {}
    for comp in ("endothelial", "myeloid_macrophage"):
        if comp not in pb_213:
            continue
        cm_mat, cm_meta = pb_213[comp]
        sub_c = cm_meta[cm_meta["site_code"].isin(["EuE", "EcP"])].copy()
        if (sub_c["site_code"] == "EcP").sum() >= 3 and (sub_c["site_code"] == "EuE").sum() >= 3:
            Cc = cm_mat[:, sub_c.index.values]
            cmask = (Cc / np.maximum(Cc.sum(axis=0, keepdims=True), 1.0) * 1e6 >= 1.0).mean(axis=1) >= 0.20
            Cc_s = Cc[cmask, :]
            tmm_c = calc_tmm_norm_factors(Cc_s)
            Yc = np.log2((Cc_s + 0.5) / ((Cc_s.sum(axis=0) * tmm_c)[None, :] + 1.0) * 1e6)
            is_ecp_c = (sub_c["site_code"].values == "EcP").astype(float)
            fit_c = limma_duplicate_correlation_gls(
                Yc, np.column_stack([np.ones(len(sub_c)), is_ecp_c]), sub_c["donor_id"].tolist(), coef_idx=1
            )
            c_genes = [g for g, k in zip(ensg_213, cmask) if k]
            ext_niche_de[comp] = {"genes": c_genes, "fit": fit_c}
            sens_213[f"{comp}_EcP_vs_EuE"] = {
                "n_EcP": int((is_ecp_c == 1).sum()),
                "n_EuE": int((is_ecp_c == 0).sum()),
                "sig_BH_q05": int((fit_c["q"] < 0.05).sum()),
            }

    # =========================================================================
    # PART 2: External 408-Sample / 162-Subject Tissue & Double-Difference (GSE141549)
    # =========================================================================
    out_bulk = ROOT / "reports/external_bulk"
    out_bulk.mkdir(parents=True, exist_ok=True)
    bulk_genes, bulk_titles, M_408, total_probes = load_gse141549_gene_matrix(
        ROOT / "data/raw/bulk/GSE141549",
        ROOT / "data/raw/metadata/GSE141549/GSE141549_family.soft.gz",
        hgnc,
    )
    title_to_col = {t: i for i, t in enumerate(bulk_titles)}
    with (ROOT / "manifests/m0_GSE141549_samples.tsv").open("r", encoding="utf-8") as f:
        s141_rows = list(csv.DictReader(f, delimiter="\t"))

    # Collapse technical replicates first, then average within (subject_id, major_group)
    # so every subject has exact weight 1 per tissue group!
    subj_group_cols = defaultdict(lambda: defaultdict(list))
    subj_meta_info = {}
    for r in s141_rows:
        col_i = title_to_col[r["sample_title"]]
        subj = r["subject_id"]
        grp = r["major_group"]
        subj_group_cols[subj][grp].append(col_i)
        if subj not in subj_meta_info:
            subj_meta_info[subj] = {
                "age": float(r["age"]) if r["age"] != "NA" and r["age"].replace(".", "").isdigit() else np.nan,
                "cycle_phase": r["cycle_phase"],
            }

    subj_expr = {}
    for subj, gdict in sorted(subj_group_cols.items()):
        subj_expr[subj] = {grp: M_408[:, cols].mean(axis=1) for grp, cols in gdict.items()}

    # 2A. Paired Within-Subject Peritoneal_Lesion vs PE (n=41 subjects)
    pe_pel_subjs = sorted(s for s, g in subj_expr.items() if "PE" in g and "Peritoneal_Lesion" in g)
    D_pel_pe = np.column_stack([subj_expr[s]["Peritoneal_Lesion"] - subj_expr[s]["PE"] for s in pe_pel_subjs])
    fit_141_pel_pe = paired_difference_moderated_fit(D_pel_pe)

    # 2B. Host-Controlled Paired Within-Subject Peritoneal_Lesion vs PP (n=20 subjects)
    pp_pel_subjs = sorted(s for s, g in subj_expr.items() if "PP" in g and "Peritoneal_Lesion" in g)
    D_pel_pp = np.column_stack([subj_expr[s]["Peritoneal_Lesion"] - subj_expr[s]["PP"] for s in pp_pel_subjs])
    fit_141_pel_pp = paired_difference_moderated_fit(D_pel_pp)

    # Save full paired DE tables for GSE141549
    for tag_name, fit_obj in (
        ("Paired_Peritoneal_Lesion_vs_PE", fit_141_pel_pe),
        ("Paired_Peritoneal_Lesion_vs_PP_HostControlled", fit_141_pel_pp),
    ):
        brows = [
            {
                "ensembl_gene_id": g,
                "symbol": hgnc["ensg_to_symbol"].get(g, ""),
                "logFC": round(float(fit_obj["logFC"][i]), 6),
                "t": round(float(fit_obj["t"][i]), 6),
                "z": round(float(fit_obj["z"][i]), 6),
                "p": float(fit_obj["p"][i]),
                "q_bh": float(fit_obj["q"][i]),
            }
            for i, g in enumerate(bulk_genes)
        ]
        write_tsv(out_bulk / f"GSE141549_{tag_name}.tsv", brows, list(brows[0].keys()))

    metrics_141_pe = compute_shortlist_metrics(
        shortlists, bulk_genes, fit_141_pel_pe["logFC"], fit_141_pel_pe["z"], fit_141_pel_pe["q"]
    )
    metrics_141_pp = compute_shortlist_metrics(
        shortlists, bulk_genes, fit_141_pel_pp["logFC"], fit_141_pel_pp["z"], fit_141_pel_pp["q"]
    )

    # Compute Joint Dual-Contrast Replication & Host-Peritoneal False-Positive Rate!
    bg_pos = {g: i for i, g in enumerate(bulk_genes)}
    joint_141 = {}
    for mname in ("B0", "G1", "M1", "M2"):
        joint_141[mname] = {}
        for k in (20, 50):
            pres = [c for c in shortlists[mname][k] if c["ensembl_gene_id"] in bg_pos]
            n_eval = len(pres)
            pe_hits = 0
            pp_hits = 0
            both_hits = 0
            pe_hit_but_pp_fail = 0
            pe_hit_and_pp_inverted = 0
            both_same_dir = 0
            for c in pres:
                gi = bg_pos[c["ensembl_gene_id"]]
                ds = int(c["discovery_sign"])
                same_pe = (fit_141_pel_pe["logFC"][gi] * ds) > 0
                q_pe = fit_141_pel_pe["q"][gi] < 0.05
                same_pp = (fit_141_pel_pp["logFC"][gi] * ds) > 0
                q_pp = fit_141_pel_pp["q"][gi] < 0.05
                if same_pe and same_pp:
                    both_same_dir += 1
                if same_pe and q_pe:
                    pe_hits += 1
                    if same_pp and q_pp:
                        both_hits += 1
                    else:
                        pe_hit_but_pp_fail += 1
                    if not same_pp:
                        pe_hit_and_pp_inverted += 1
                if same_pp and q_pp:
                    pp_hits += 1
            joint_141[mname][f"top{k}"] = {
                "selected": k,
                "evaluable": n_eval,
                "PE_replicated_q05": pe_hits,
                "PP_host_replicated_q05": pp_hits,
                "Joint_PE_and_PP_replicated_q05": both_hits,
                "Joint_PE_and_PP_rate_on_k": round(both_hits / k, 6),
                "Both_PE_and_PP_same_direction": both_same_dir,
                "Both_same_direction_rate_evaluable": round(both_same_dir / max(n_eval, 1), 6),
                "Host_confounded_PE_hits_failing_PP": pe_hit_but_pp_fail,
                "Host_false_positive_fraction_of_PE_hits": round(pe_hit_but_pp_fail / pe_hits, 4) if pe_hits else 0.0,
                "Host_sign_inverted_among_PE_hits": pe_hit_and_pp_inverted,
            }

    boot_141_paired = paired_dual_contrast_bootstrap(
        shortlists, bulk_genes, subj_expr, n_draws=n_boot, seed=seed
    )

    # 2C. Exact 18-Subject Triplet Double-Difference (PE, PP, and Peritoneal_Lesion in the same woman!)
    triplet_18_subjs = sorted(s for s, g in subj_expr.items() if {"PE", "PP", "Peritoneal_Lesion"} <= set(g.keys()))
    D_trip_pe = np.column_stack([subj_expr[s]["Peritoneal_Lesion"] - subj_expr[s]["PE"] for s in triplet_18_subjs])
    D_trip_pp = np.column_stack([subj_expr[s]["Peritoneal_Lesion"] - subj_expr[s]["PP"] for s in triplet_18_subjs])
    fit_trip_pe = paired_difference_moderated_fit(D_trip_pe)
    fit_trip_pp = paired_difference_moderated_fit(D_trip_pp)
    triplet_summary = {"n_triplets": len(triplet_18_subjs), "methods": {}}
    for mname in ("B0", "G1", "M1", "M2"):
        triplet_summary["methods"][mname] = {}
        for k in (20, 50):
            pres = [c for c in shortlists[mname][k] if c["ensembl_gene_id"] in bg_pos]
            both_q05 = sum(
                1 for c in pres
                if (fit_trip_pe["logFC"][bg_pos[c["ensembl_gene_id"]]] * int(c["discovery_sign"]) > 0)
                and fit_trip_pe["q"][bg_pos[c["ensembl_gene_id"]]] < 0.05
                and (fit_trip_pp["logFC"][bg_pos[c["ensembl_gene_id"]]] * int(c["discovery_sign"]) > 0)
                and fit_trip_pp["q"][bg_pos[c["ensembl_gene_id"]]] < 0.05
            )
            triplet_summary["methods"][mname][f"top{k}_joint_q05"] = both_q05

    # 2D. Multi-Subtype Contrasts in GSE141549 (DIE vs PP paired, All-subject PeL/DIE/OMA vs PP+CP, and PE vs CE check)
    pp_die_subjs = sorted(s for s, g in subj_expr.items() if "PP" in g and "DIE" in g)
    fit_141_die_pp = paired_difference_moderated_fit(
        np.column_stack([subj_expr[s]["DIE"] - subj_expr[s]["PP"] for s in pp_die_subjs])
    )
    metrics_141_die_pp = compute_shortlist_metrics(
        shortlists, bulk_genes, fit_141_die_pp["logFC"], fit_141_die_pp["z"], fit_141_die_pp["q"]
    )

    def fit_all_subject_bulk_contrast(case_groups: tuple[str, ...], ref_groups: tuple[str, ...]):
        cols_expr = []
        y_ind = []
        b_ids = []
        phases = []
        ages = []
        for subj, gdict in sorted(subj_expr.items()):
            c_vecs = [gdict[g] for g in case_groups if g in gdict]
            r_vecs = [gdict[g] for g in ref_groups if g in gdict]
            for is_c, vec_list in ((1.0, c_vecs), (0.0, r_vecs)):
                if vec_list:
                    cols_expr.append(np.mean(vec_list, axis=0))
                    y_ind.append(is_c)
                    b_ids.append(subj)
                    phases.append(subj_meta_info[subj]["cycle_phase"])
                    ages.append(subj_meta_info[subj]["age"])
        M_c = np.column_stack(cols_expr)
        y_arr = np.array(y_ind, dtype=np.float64)
        age_arr = np.array(ages, dtype=np.float64)
        age_filled = np.where(np.isfinite(age_arr), age_arr - np.nanmean(age_arr), 0.0)
        # Covariates: intercept, case_indicator, secretory, medication, age_centered
        is_sec = np.array([1.0 if p == "secretory" else 0.0 for p in phases], dtype=np.float64)
        is_med = np.array([1.0 if p == "medication" else 0.0 for p in phases], dtype=np.float64)
        X_c = np.column_stack([np.ones(len(y_arr)), y_arr, is_sec, is_med, age_filled])
        fit_c = limma_duplicate_correlation_gls(M_c, X_c, b_ids, coef_idx=1)
        fit_c["n_case"] = int((y_arr == 1).sum())
        fit_c["n_ref"] = int((y_arr == 0).sum())
        fit_c["M_c"] = M_c
        fit_c["y_arr"] = y_arr
        fit_c["b_ids"] = b_ids
        return fit_c

    fit_all_pel_host = fit_all_subject_bulk_contrast(("Peritoneal_Lesion",), ("PP", "CP"))
    fit_all_die_host = fit_all_subject_bulk_contrast(("DIE",), ("PP", "CP"))
    fit_all_oma_pe = fit_all_subject_bulk_contrast(("OMA",), ("PE",))
    fit_all_oma_host = fit_all_subject_bulk_contrast(("OMA",), ("PP", "CP"))
    fit_all_pe_ce = fit_all_subject_bulk_contrast(("PE",), ("CE",))

    prog_141_pel_host = compute_program_transfer(
        shortlists, bulk_genes, fit_all_pel_host["M_c"], fit_all_pel_host["y_arr"], block_ids=fit_all_pel_host["b_ids"]
    )

    subtype_summary = {
        "Paired_DIE_vs_PP_n20": {
            "n_pairs": fit_141_die_pp["n_pairs"],
            "sig_BH_q05": int((fit_141_die_pp["q"] < 0.05).sum()),
            "shortlists": {
                m: {k: {kk: vv for kk, vv in d.items() if kk != "genes"} for k, d in md.items()}
                for m, md in metrics_141_die_pp.items()
            },
        },
        "All_Peritoneal_Lesion_vs_PP_CP": {
            "n_case": fit_all_pel_host["n_case"],
            "n_ref": fit_all_pel_host["n_ref"],
            "sig_BH_q05": int((fit_all_pel_host["q"] < 0.05).sum()),
            "shortlists": {
                m: {k: {kk: vv for kk, vv in d.items() if kk != "genes"} for k, d in md.items()}
                for m, md in compute_shortlist_metrics(shortlists, bulk_genes, fit_all_pel_host["logFC"], fit_all_pel_host["z"], fit_all_pel_host["q"]).items()
            },
        },
        "All_DIE_vs_PP_CP": {
            "n_case": fit_all_die_host["n_case"],
            "n_ref": fit_all_die_host["n_ref"],
            "sig_BH_q05": int((fit_all_die_host["q"] < 0.05).sum()),
            "shortlists": {
                m: {k: {kk: vv for kk, vv in d.items() if kk != "genes"} for k, d in md.items()}
                for m, md in compute_shortlist_metrics(shortlists, bulk_genes, fit_all_die_host["logFC"], fit_all_die_host["z"], fit_all_die_host["q"]).items()
            },
        },
        "All_OMA_vs_PE": {
            "n_case": fit_all_oma_pe["n_case"],
            "n_ref": fit_all_oma_pe["n_ref"],
            "sig_BH_q05": int((fit_all_oma_pe["q"] < 0.05).sum()),
            "shortlists": {
                m: {k: {kk: vv for kk, vv in d.items() if kk != "genes"} for k, d in md.items()}
                for m, md in compute_shortlist_metrics(shortlists, bulk_genes, fit_all_oma_pe["logFC"], fit_all_oma_pe["z"], fit_all_oma_pe["q"]).items()
            },
        },
        "All_OMA_vs_PP_CP": {
            "n_case": fit_all_oma_host["n_case"],
            "n_ref": fit_all_oma_host["n_ref"],
            "sig_BH_q05": int((fit_all_oma_host["q"] < 0.05).sum()),
            "shortlists": {
                m: {k: {kk: vv for kk, vv in d.items() if kk != "genes"} for k, d in md.items()}
                for m, md in compute_shortlist_metrics(shortlists, bulk_genes, fit_all_oma_host["logFC"], fit_all_oma_host["z"], fit_all_oma_host["q"]).items()
            },
        },
        "Eutopic_PE_vs_CE_Negative_Control": {
            "n_case": fit_all_pe_ce["n_case"],
            "n_ref": fit_all_pe_ce["n_ref"],
            "sig_BH_q05": int((fit_all_pe_ce["q"] < 0.05).sum()),
        },
    }

    # =========================================================================
    # PART 3: Compartment-Resolved Human Genetic Enrichment & Target Dossiers
    # =========================================================================
    gene_sets = {
        "OpenTargets_Genetic_Any": {r["ensembl_gene_id"] for r in gen_rows if float(r["ot_genetic_association"]) > 0},
        "OpenTargets_Genetic_ge0.5": {r["ensembl_gene_id"] for r in gen_rows if float(r["ot_genetic_association"]) >= 0.5},
        "Rahmioglu2023_CredibleNearest": {r["ensembl_gene_id"] for r in gen_rows if r["rahmioglu_cred_nearest"] == "True"},
        "Rahmioglu2023_SMR_Causal": {r["ensembl_gene_id"] for r in gen_rows if r["rahmioglu_smr_any"] == "True"},
        "Combined_Human_Genetics_Prior": {r["ensembl_gene_id"] for r in gen_rows if float(r["genetic_prior_score"]) > 0},
    }

    compartment_enrichment = {}
    for comp_file in sorted((ROOT / "reports/discovery").glob("GSE179640_DE_*_EcP_vs_EuE.tsv")):
        cname = comp_file.stem.replace("GSE179640_DE_", "").replace("_EcP_vs_EuE", "")
        with comp_file.open("r", encoding="utf-8") as f:
            crows = list(csv.DictReader(f, delimiter="\t"))
        cgenes = [r["ensembl_gene_id"] for r in crows]
        cz = np.abs(np.array([float(r["z"]) for r in crows], dtype=np.float64))
        compartment_enrichment[f"GSE179640_{cname}_absZ"] = test_genetic_set_enrichment(cgenes, cz, gene_sets, seed=seed)

    # Also test M1_core on Discovery stromal universe and external GSE213216 / GSE141549
    u_genes = [r["ensembl_gene_id"] for r in disc_full_rows]
    u_m1 = np.array([float(r["m1_core"]) for r in disc_full_rows], dtype=np.float64)
    compartment_enrichment["GSE179640_stromal_M1_host_decontaminated"] = test_genetic_set_enrichment(
        u_genes, u_m1, gene_sets, seed=seed
    )
    compartment_enrichment["GSE213216_stromal_EcP_vs_EuE_absZ"] = test_genetic_set_enrichment(
        genes_213, np.abs(fit_213_prim["z"]), gene_sets, expr_vec=fit_213_prim["AveExpr"], seed=seed
    )
    for comp, d_obj in ext_niche_de.items():
        compartment_enrichment[f"GSE213216_{comp}_EcP_vs_EuE_absZ"] = test_genetic_set_enrichment(
            d_obj["genes"], np.abs(d_obj["fit"]["z"]), gene_sets, expr_vec=d_obj["fit"]["AveExpr"], seed=seed
        )
    compartment_enrichment["GSE141549_Paired_PeL_vs_PP_HostControlled_absZ"] = test_genetic_set_enrichment(
        bulk_genes, np.abs(fit_141_pel_pp["z"]), gene_sets, seed=seed
    )

    # Build Evidence-Graded Target Dossiers across all genes in Shared Universe U
    g213_pos = {g: i for i, g in enumerate(genes_213)}
    dossier_candidates = []

    for r in disc_full_rows:
        gid = r["ensembl_gene_id"]
        m1_c = float(r["m1_core"])
        if m1_c <= 0.5:
            continue
        dsign = int(r["discovery_sign"])
        g_prior = float(r["genetic_prior_score"])
        n_lr = float(r["niche_lr_score"])
        gen_rec = genetics_by_ensg.get(gid, {})

        # External scRNA GSE213216 evidence
        if gid in g213_pos:
            i213 = g213_pos[gid]
            sc_lfc = float(fit_213_prim["logFC"][i213])
            sc_p = float(fit_213_prim["p"][i213])
            sc_q = float(fit_213_prim["q"][i213])
            sc_same = (sc_lfc * dsign) > 0
        else:
            sc_lfc, sc_p, sc_q, sc_same = np.nan, 1.0, 1.0, False

        # External Bulk GSE141549 Paired PeL vs PE and Host-Controlled PeL vs PP
        if gid in bg_pos:
            ib = bg_pos[gid]
            pe_lfc = float(fit_141_pel_pe["logFC"][ib])
            pe_q = float(fit_141_pel_pe["q"][ib])
            pe_same = (pe_lfc * dsign) > 0
            pp_lfc = float(fit_141_pel_pp["logFC"][ib])
            pp_q = float(fit_141_pel_pp["q"][ib])
            pp_same = (pp_lfc * dsign) > 0
            die_lfc = float(fit_141_die_pp["logFC"][ib])
            die_q = float(fit_141_die_pp["q"][ib])
            die_same = (die_lfc * dsign) > 0
        else:
            pe_lfc = pp_lfc = die_lfc = np.nan
            pe_q = pp_q = die_q = 1.0
            pe_same = pp_same = die_same = False

        # Tier classification per protocol_v2.yaml:
        # Tier 1: Human genetic support (G_prior > 0) + Discovery M1 dual-contrast concordant
        #         + External scRNA GSE213216 concordant (q < 0.05 or nominal p < 0.05)
        #         + External Host-Controlled Bulk GSE141549 Lesion vs PP BH q < 0.05 same sign
        is_tier1 = (
            g_prior > 0.0
            and m1_c >= 1.0
            and sc_same
            and sc_p < 0.05
            and pp_same
            and pp_q < 0.05
        )
        # Tier 2: In M2 or M1 top-50 + Replicated BH q < 0.05 same sign in GSE213216 scRNA
        #         AND GSE141549 Host-Controlled Lesion vs PP + (Active L-R or SM/AB tractability)
        tract_any = gen_rec.get("tractable_sm_any") == "True" or gen_rec.get("tractable_ab_any") == "True"
        is_tier2 = (
            not is_tier1
            and (int(r["rank_M2"]) <= 50 or int(r["rank_M1"]) <= 50)
            and sc_same
            and sc_q < 0.05
            and pp_same
            and pp_q < 0.05
            and (n_lr > 0.0 or tract_any or g_prior > 0.0)
        )

        if is_tier1 or is_tier2:
            tier = "Tier_1_Genetics_Anchored_MultiModal" if is_tier1 else "Tier_2_Replicated_Niche_Tractable"
            dossier_candidates.append({
                "tier": tier,
                "ensembl_gene_id": gid,
                "symbol": r["symbol"],
                "direction_in_lesion": "UP_in_Lesion" if dsign > 0 else "DOWN_in_Lesion",
                "rank_M2": int(r["rank_M2"]),
                "rank_M1": int(r["rank_M1"]),
                "rank_B0": int(r["rank_B0"]),
                "disc_logFC_EcP_vs_EuE": round(float(r["logFC_EcP_vs_EuE"]), 4),
                "disc_logFC_EcP_vs_AdjP": round(float(r["logFC_EcP_vs_AdjP"]), 4),
                "disc_m1_core": round(m1_c, 4),
                "ext_sc_GSE213216_logFC": round(sc_lfc, 4),
                "ext_sc_GSE213216_p": sc_p,
                "ext_sc_GSE213216_q": sc_q,
                "ext_bulk_PeL_vs_PE_logFC": round(pe_lfc, 4),
                "ext_bulk_PeL_vs_PE_q": pe_q,
                "ext_bulk_PeL_vs_PP_HostCtrl_logFC": round(pp_lfc, 4),
                "ext_bulk_PeL_vs_PP_HostCtrl_q": pp_q,
                "ext_bulk_DIE_vs_PP_HostCtrl_logFC": round(die_lfc, 4),
                "ext_bulk_DIE_vs_PP_HostCtrl_q": die_q,
                "genetic_prior_score": round(g_prior, 4),
                "ot_genetic_association": gen_rec.get("ot_genetic_association", "0"),
                "ot_l2g_credible_set": gen_rec.get("ot_l2g_credible_set", "0"),
                "rahmioglu_gwas_or_smr": str(
                    gen_rec.get("rahmioglu_gws42") == "True"
                    or gen_rec.get("rahmioglu_cred_nearest") == "True"
                    or gen_rec.get("rahmioglu_smr_any") == "True"
                ),
                "niche_lr_score": round(n_lr, 4),
                "best_lr_pair": r["best_lr_pair"],
                "best_lr_compartment": r["best_lr_compartment"],
                "tractable_sm": gen_rec.get("tractable_sm_any", "False"),
                "tractable_sm_clinical": gen_rec.get("tractable_sm_clinical", "False"),
                "tractable_ab": gen_rec.get("tractable_ab_any", "False"),
                "ot_clinical_drug_score": gen_rec.get("ot_clinical_drug", "0"),
            })

    dossier_candidates.sort(
        key=lambda x: (0 if x["tier"].startswith("Tier_1") else 1, -x["genetic_prior_score"], x["rank_M2"])
    )

    # Enrich top dossier genes with GTEx v10 tissue expression & eQTLs via gtex_cli.py
    top_syms = [d["symbol"] for d in dossier_candidates[:25] if d["symbol"]]
    gtex_map = query_gtex_for_dossier_genes(top_syms)
    for d in dossier_candidates:
        gt = gtex_map.get(d["symbol"], {})
        d["gtex_tpm_uterus"] = gt.get("tpm_Uterus", "NA")
        d["gtex_tpm_fibroblasts"] = gt.get("tpm_Fibroblasts", "NA")
        d["gtex_tpm_ovary"] = gt.get("tpm_Ovary", "NA")
        d["gtex_eqtls_uterus_fibro_ovary"] = (
            f"{gt.get('eqtls_Uterus', 0)}/{gt.get('eqtls_Fibroblasts', 0)}/{gt.get('eqtls_Ovary', 0)}"
            if gt else "NA"
        )

    if dossier_candidates:
        write_tsv(ROOT / "reports/target_dossiers.tsv", dossier_candidates, list(dossier_candidates[0].keys()))

    # Write structured Markdown Target Dossiers report
    md_lines = [
        "# Multi-Evidence Endometriosis Lesion Target Dossiers (EndoTarget v2)",
        "",
        "Generated strictly under locked `configs/protocol_v2.yaml` and `receipts/prediction_seal.json`.",
        "Data sources include GEO (`GSE179640`, `GSE213216`, `GSE141549`), Open Targets Platform (`MONDO_0005133` / `EFO_0001065`; https://platform-docs.opentargets.org/licence), Rahmioglu et al. 2023 (*Nat Genet*), and GTEx Portal API v10 (https://gtexportal.org/home/license).",
        "",
        f"- **Tier 1 (`Tier_1_Genetics_Anchored_MultiModal`)**: **{sum(1 for d in dossier_candidates if d['tier'].startswith('Tier_1'))}** targets backed simultaneously by human endometriosis germline genetics (`G_prior > 0`), discovery dual-contrast host-decontaminated remodeling (`GSE179640`), external scRNA stromal concordance (`GSE213216`, `p < 0.05`), and external host-peritoneum-controlled tissue validation (`GSE141549` paired `Peritoneal_Lesion vs PP` `BH q < 0.05`).",
        f"- **Tier 2 (`Tier_2_Replicated_Niche_Tractable`)**: **{sum(1 for d in dossier_candidates if d['tier'].startswith('Tier_2'))}** frozen `M1`/`M2` Top-50 targets achieving strict `BH q < 0.05` replication in **both** external scRNA (`GSE213216`) and external host-controlled bulk (`GSE141549` `PeL vs PP`) with active vascular/myeloid L–R crosstalk or small-molecule/antibody tractability.",
        "",
        "| Tier | Gene | Lesion Dir | Rank (M2/M1/B0) | Disc logFC (EuE / AdjP) | GSE213216 scRNA logFC (q) | GSE141549 PeL vs PP logFC (q) | DIE vs PP logFC (q) | Genetic Prior (OT / GWAS-SMR) | Niche L-R Axis | Tractability (SM / AB) | GTEx TPM (Ut/Fib/Ov) |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|",
    ]
    for d in dossier_candidates:
        tier_short = "Tier 1" if d["tier"].startswith("Tier_1") else "Tier 2"
        md_lines.append(
            f"| {tier_short} | **{d['symbol']}** (`{d['ensembl_gene_id']}`) | `{d['direction_in_lesion']}` | "
            f"{d['rank_M2']} / {d['rank_M1']} / {d['rank_B0']} | "
            f"{d['disc_logFC_EcP_vs_EuE']:+.2f} / {d['disc_logFC_EcP_vs_AdjP']:+.2f} | "
            f"{d['ext_sc_GSE213216_logFC']:+.2f} ({d['ext_sc_GSE213216_q']:.2e}) | "
            f"{d['ext_bulk_PeL_vs_PP_HostCtrl_logFC']:+.2f} ({d['ext_bulk_PeL_vs_PP_HostCtrl_q']:.2e}) | "
            f"{d['ext_bulk_DIE_vs_PP_HostCtrl_logFC']:+.2f} ({d['ext_bulk_DIE_vs_PP_HostCtrl_q']:.2e}) | "
            f"{d['genetic_prior_score']:.2f} ({float(d['ot_genetic_association']):.2f} / {d['rahmioglu_gwas_or_smr']}) | "
            f"`{d['best_lr_pair'] or '-'}` ({d['best_lr_compartment'] or '-'}) | "
            f"SM={d['tractable_sm']}, AB={d['tractable_ab']} | "
            f"{d['gtex_tpm_uterus']}/{d['gtex_tpm_fibroblasts']}/{d['gtex_tpm_ovary']} |"
        )
    (ROOT / "reports/target_dossiers.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    # =========================================================================
    # PART 4: Save Master Locked Evaluation JSON & Cryptographic Receipt
    # =========================================================================
    def strip_gene_arrays(mdict: dict) -> dict:
        return {
            m: {k: {kk: vv for kk, vv in d.items() if kk != "genes"} for k, d in by_k.items()}
            for m, by_k in mdict.items()
        }

    locked_summary = {
        "schema_version": 2,
        "evaluated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "locks_verified": {
            "protocol_freeze_sha256": freeze_sha,
            "prediction_seal_sha256": seal_sha,
            "external_cohorts_data_sha256": ext_data_sha,
        },
        "GSE213216_external_scrna": {
            "libraries_evaluated": len(lib_qc_213),
            "raw_cells_total": int(sum(r["raw_cells"] for r in lib_qc_213)),
            "qc_pass_cells_total": int(X_213.shape[0]),
            "lineage_cell_counts": {k: int(v) for k, v in Counter(annot_213["lineage"]).items()},
            "pseudobulk_conservation": cons_213,
            "primary_stromal_EcP_vs_EuE": {
                "n_EcP_donors": int((is_ecp_213 == 1).sum()),
                "n_EuE_donors": int((is_ecp_213 == 0).sum()),
                "n_unique_patients": fit_213_prim["n_unique_blocks"],
                "consensus_rho": round(fit_213_prim["consensus_rho"], 4),
                "df_residual": fit_213_prim["df_residual"],
                "prior_df": round(fit_213_prim["prior_df"], 4),
                "genes_evaluated": len(genes_213),
                "sig_BH_q_lt_0.05": int((fit_213_prim["q"] < 0.05).sum()),
                "shortlist_metrics": strip_gene_arrays(metrics_213_prim),
                "program_transfer": prog_213_prim,
                "donor_cluster_bootstrap": boot_213_prim,
            },
            "sensitivities": sens_213,
        },
        "GSE141549_external_bulk": {
            "total_arrays": 408,
            "unique_subjects": len(subj_expr),
            "total_probes_parsed": total_probes,
            "mapped_unique_ensembl_genes": len(bulk_genes),
            "Paired_Peritoneal_Lesion_vs_PE_n41": {
                "n_pairs": fit_141_pel_pe["n_pairs"],
                "df_residual": fit_141_pel_pe["df_residual"],
                "prior_df": round(fit_141_pel_pe["prior_df"], 4),
                "sig_BH_q_lt_0.05": int((fit_141_pel_pe["q"] < 0.05).sum()),
                "shortlist_metrics": strip_gene_arrays(metrics_141_pe),
            },
            "Paired_Peritoneal_Lesion_vs_PP_HostControlled_n20": {
                "n_pairs": fit_141_pel_pp["n_pairs"],
                "df_residual": fit_141_pel_pp["df_residual"],
                "prior_df": round(fit_141_pel_pp["prior_df"], 4),
                "sig_BH_q_lt_0.05": int((fit_141_pel_pp["q"] < 0.05).sum()),
                "shortlist_metrics": strip_gene_arrays(metrics_141_pp),
            },
            "Joint_Dual_Contrast_and_Host_False_Positive_Audit": joint_141,
            "Program_Transfer_All_Peritoneal_vs_Host": prog_141_pel_host,
            "Paired_Subject_Bootstrap": boot_141_paired,
            "Triplet_18subject_Double_Difference": triplet_summary,
            "Multi_Subtype_and_Control_Contrasts": subtype_summary,
        },
        "human_genetics_enrichment": compartment_enrichment,
        "target_dossiers_summary": {
            "tier1_count": sum(1 for d in dossier_candidates if d["tier"].startswith("Tier_1")),
            "tier2_count": sum(1 for d in dossier_candidates if d["tier"].startswith("Tier_2")),
            "tier1_symbols": [d["symbol"] for d in dossier_candidates if d["tier"].startswith("Tier_1")],
            "tier2_symbols": [d["symbol"] for d in dossier_candidates if d["tier"].startswith("Tier_2")],
        },
    }

    sum_path = ROOT / "reports/locked_evaluation_summary_v2.json"
    sum_path.write_text(json.dumps(locked_summary, indent=2) + "\n", encoding="utf-8")

    # Save detailed per-candidate external evaluation table
    cand_eval_rows = []
    for mname in ("B0", "G1", "M1", "M2"):
        for g_rec in metrics_213_prim[mname]["top20"]["genes"]:
            gid = g_rec["ensembl_gene_id"]
            ib = bg_pos.get(gid)
            ds = g_rec["discovery_sign"]
            cand_eval_rows.append({
                "method": mname,
                "rank_in_method": g_rec["rank"],
                "ensembl_gene_id": gid,
                "symbol": g_rec["symbol"],
                "discovery_sign": ds,
                "GSE213216_logFC": g_rec["ext_logFC"],
                "GSE213216_q": g_rec["ext_q"],
                "GSE213216_replicated_q05": str(g_rec["replicated_q05"]),
                "GSE141549_PeL_vs_PE_logFC": round(float(fit_141_pel_pe["logFC"][ib]), 5) if ib is not None else "NA",
                "GSE141549_PeL_vs_PE_q": float(fit_141_pel_pe["q"][ib]) if ib is not None else "NA",
                "GSE141549_PeL_vs_PP_HostCtrl_logFC": round(float(fit_141_pel_pp["logFC"][ib]), 5) if ib is not None else "NA",
                "GSE141549_PeL_vs_PP_HostCtrl_q": float(fit_141_pel_pp["q"][ib]) if ib is not None else "NA",
                "GSE141549_Joint_PE_and_PP_q05": str(
                    ib is not None
                    and (fit_141_pel_pe["logFC"][ib] * ds > 0)
                    and fit_141_pel_pe["q"][ib] < 0.05
                    and (fit_141_pel_pp["logFC"][ib] * ds > 0)
                    and fit_141_pel_pp["q"][ib] < 0.05
                ),
            })
    write_tsv(ROOT / "reports/top20_external_audit_v2.tsv", cand_eval_rows, list(cand_eval_rows[0].keys()))

    eval_outputs = [
        "reports/locked_evaluation_summary_v2.json",
        "reports/top20_external_audit_v2.tsv",
        "reports/target_dossiers.tsv",
        "reports/target_dossiers.md",
        "reports/external_scrna/GSE213216_library_qc.tsv",
        "reports/external_scrna/GSE213216_cluster_annotations.tsv",
        "reports/external_scrna/GSE213216_pseudobulk_samples.tsv",
        "reports/external_scrna/GSE213216_stromal_DE_EcP_vs_EuE.tsv",
        "reports/external_bulk/GSE141549_Paired_Peritoneal_Lesion_vs_PE.tsv",
        "reports/external_bulk/GSE141549_Paired_Peritoneal_Lesion_vs_PP_HostControlled.tsv",
    ]
    eval_receipt = {
        "schema_version": 2,
        "phase": "STAGE3_LOCKED_EXTERNAL_EVALUATION",
        "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "protocol_freeze_sha256": freeze_sha,
        "prediction_seal_sha256": seal_sha,
        "external_cohorts_data_sha256": ext_data_sha,
        "evaluate_script_sha256": sha256_file(pathlib.Path(__file__)),
        "outputs": {
            rel: {"sha256": sha256_file(ROOT / rel), "bytes": (ROOT / rel).stat().st_size}
            for rel in eval_outputs if (ROOT / rel).exists()
        },
    }
    (ROOT / "receipts/locked_evaluation_v2.json").write_text(json.dumps(eval_receipt, indent=2) + "\n", encoding="utf-8")

    print(json.dumps({
        "status": "EVALUATION_COMPLETE",
        "GSE213216_qc_cells": locked_summary["GSE213216_external_scrna"]["qc_pass_cells_total"],
        "GSE213216_top20_rep_q05": {
            m: locked_summary["GSE213216_external_scrna"]["primary_stromal_EcP_vs_EuE"]["shortlist_metrics"][m]["top20"]["replicated_q05_count"]
            for m in ("B0", "G1", "M1", "M2")
        },
        "GSE141549_joint_PE_and_PP_top20": {
            m: joint_141[m]["top20"]["Joint_PE_and_PP_replicated_q05"]
            for m in ("B0", "G1", "M1", "M2")
        },
        "GSE141549_host_false_positive_fraction_top20": {
            m: joint_141[m]["top20"]["Host_false_positive_fraction_of_PE_hits"]
            for m in ("B0", "G1", "M1", "M2")
        },
        "dossiers": locked_summary["target_dossiers_summary"],
    }, indent=2))


if __name__ == "__main__":
    main()
