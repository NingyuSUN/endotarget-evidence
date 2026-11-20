#!/usr/bin/env python3
"""Stage 2: Discovery single-cell processing (GSE179640), de novo canonical lineage
pseudobulk, dual-contrast host-decontaminated modeling (EcP vs EuE & EcP vs AdjP),
niche L-R scoring, 4-method candidate freezing (B0, G1, M1, M2), and prediction_seal.json.

STRICT ZERO-LEAKAGE GUARANTEE:
- Verifies receipts/protocol_freeze.json and receipts/discovery_GSE179640_data.json before running.
- Asserts zero external cohort expression files exist in data/raw/scrna/GSE213216 or data/raw/bulk.
"""
from __future__ import annotations

import csv
import datetime as dt
import gzip
import json
import math
import os
import pathlib
import sys
from collections import Counter, defaultdict

os.environ["OMP_NUM_THREADS"] = "8"
os.environ["OPENBLAS_NUM_THREADS"] = "8"
os.environ["MKL_NUM_THREADS"] = "8"
os.environ["NUMEXPR_NUM_THREADS"] = "8"

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from endotarget_lib import (  # noqa: E402
    aggregate_pseudobulk,
    assign_lineages_de_novo,
    limma_duplicate_correlation_gls,
    limma_voom_fit,
    load_10x_h5_canonical,
    rwr_propagate,
)
from freeze_protocol import sha256_file, verify_protocol_freeze  # noqa: E402
from m0_audit_and_priors import load_hgnc, write_tsv  # noqa: E402


def align_sparse_matrices(mats: list[sp.csr_matrix], gene_lists: list[list[str]]) -> tuple[sp.csr_matrix, list[str]]:
    common = sorted(set().union(*gene_lists))
    if all(gl == common for gl in gene_lists):
        return sp.vstack(mats, format="csr"), common
    pos = {g: i for i, g in enumerate(common)}
    aligned = []
    for X, gl in zip(mats, gene_lists):
        cols = np.array([pos[g] for g in gl], dtype=np.int64)
        proj = sp.csr_matrix(
            (np.ones(len(gl), dtype=np.int64), (np.arange(len(gl), dtype=np.int64), cols)),
            shape=(len(gl), len(common)),
        )
        aligned.append((X @ proj).tocsr())
    return sp.vstack(aligned, format="csr"), common


def fit_paired_donor_voom(counts: np.ndarray, meta_df: pd.DataFrame, site_case: str, site_ref: str, keep_mask: np.ndarray | None = None):
    """Fit paired ~ donor + is_case limma-voom model on donors possessing both site_case and site_ref."""
    sub_meta = meta_df[meta_df["site_code"].isin([site_ref, site_case])].copy()
    by_donor = sub_meta.groupby("donor_id")["site_code"].apply(set)
    paired_donors = sorted(d for d, s in by_donor.items() if {site_ref, site_case} <= s)
    if len(paired_donors) < 3:
        raise ValueError(f"Too few paired donors ({len(paired_donors)}) for {site_case} vs {site_ref}")

    sel = sub_meta[sub_meta["donor_id"].isin(paired_donors)].sort_values(["donor_id", "site_code"])
    col_indices = sel.index.values
    C = counts[:, col_indices]
    if keep_mask is not None:
        C = C[keep_mask, :]

    # Design matrix: one indicator per donor (sum/dummy coding with first donor baseline + intercept) + case indicator
    n_samp = len(sel)
    donors = sel["donor_id"].tolist()
    u_donors = paired_donors
    X_cols = [np.ones(n_samp, dtype=np.float64)]
    for d in u_donors[1:]:
        X_cols.append(np.array([1.0 if x == d else 0.0 for x in donors], dtype=np.float64))
    is_case = np.array([1.0 if s == site_case else 0.0 for s in sel["site_code"]], dtype=np.float64)
    X_cols.append(is_case)
    X = np.column_stack(X_cols)

    fit = limma_voom_fit(C, X, coef_idx=X.shape[1] - 1)
    fit["paired_donors"] = u_donors
    fit["n_pairs"] = len(u_donors)
    fit["sample_meta"] = sel.reset_index(drop=True)
    return fit


def main():
    # 1. Cryptographic zero-leakage verification
    verify_protocol_freeze(ROOT)
    freeze_sha = sha256_file(ROOT / "receipts/protocol_freeze.json")
    cfg = yaml.safe_load((ROOT / "configs/protocol_v2.yaml").read_text(encoding="utf-8"))

    for ext_dir in (ROOT / "data/raw/scrna/GSE213216", ROOT / "data/raw/bulk"):
        if ext_dir.exists() and any(ext_dir.rglob("*")):
            raise RuntimeError(f"Zero-leakage violation: external data directory {ext_dir} exists during Stage 2!")

    data_rec_path = ROOT / "receipts/discovery_GSE179640_data.json"
    data_rec = json.loads(data_rec_path.read_text(encoding="utf-8"))
    if data_rec["protocol_freeze_sha256"] != freeze_sha:
        raise ValueError("Discovery data receipt protocol_freeze_sha256 mismatch")

    h5_dir = ROOT / "data/raw/scrna/GSE179640"
    for fname, rec in data_rec["files"].items():
        if sha256_file(h5_dir / fname) != rec["sha256"]:
            raise ValueError(f"SHA-256 mismatch on discovery h5 file {fname}")

    hgnc = load_hgnc(ROOT / "data/raw/reference/hgnc_complete_set.tsv")
    sym_to_ensg = {**hgnc["alias"], **hgnc["prev"], **hgnc["approved"]}

    with (ROOT / "manifests/canonical_lineage_markers.tsv").open("r", encoding="utf-8") as f:
        marker_rows = list(csv.DictReader(f, delimiter="\t"))
    with (ROOT / "manifests/curated_niche_lr_pairs.tsv").open("r", encoding="utf-8") as f:
        lr_rows = list(csv.DictReader(f, delimiter="\t"))
    with (ROOT / "manifests/human_genetics_prior_v2.tsv").open("r", encoding="utf-8") as f:
        genetics_by_ensg = {r["ensembl_gene_id"]: r for r in csv.DictReader(f, delimiter="\t")}

    # 2. Load 31 discovery .h5 matrices and apply locked single-cell QC
    qc_cfg = cfg["single_cell_processing"]["qc_filters"]
    mats = []
    glists = []
    obs_parts = []
    lib_qc_rows = []

    for fname in sorted(data_rec["files"].keys()):
        rec = data_rec["files"][fname]
        X_lib, ensg_lib, qc_df = load_10x_h5_canonical(h5_dir / fname, sym_to_ensg=sym_to_ensg)
        pass_qc = (
            (qc_df["n_genes_by_counts"] >= int(qc_cfg["min_genes_by_counts"]))
            & (qc_df["total_counts"] >= int(qc_cfg["min_total_counts"]))
            & (qc_df["pct_counts_mt"] <= float(qc_cfg["max_pct_counts_mt"]))
        ).values

        lib_qc_rows.append({
            "gsm": rec["gsm"],
            "donor_id": rec["donor_id"],
            "site_code": rec["site_code"],
            "filename": fname,
            "raw_cells": int(X_lib.shape[0]),
            "qc_pass_cells": int(pass_qc.sum()),
            "qc_pass_fraction": round(float(pass_qc.mean()), 4),
            "median_genes_qc": float(np.median(qc_df.loc[pass_qc, "n_genes_by_counts"])) if pass_qc.any() else 0.0,
            "median_umis_qc": float(np.median(qc_df.loc[pass_qc, "total_counts"])) if pass_qc.any() else 0.0,
        })

        sub_qc = qc_df.loc[pass_qc].copy()
        sub_qc["gsm"] = rec["gsm"]
        sub_qc["donor_id"] = rec["donor_id"]
        sub_qc["site_code"] = rec["site_code"]
        sub_qc["cell_id"] = rec["gsm"] + ":" + sub_qc["barcode"].astype(str)
        mats.append(X_lib[pass_qc])
        glists.append(ensg_lib)
        obs_parts.append(sub_qc)

    X_all, ensg_ids = align_sparse_matrices(mats, glists)
    obs_all = pd.concat(obs_parts, ignore_index=True)
    obs_all.index = obs_all["cell_id"].values

    # 3. De novo canonical lineage assignment (strictly zero HECA labels!)
    annot_obs, cluster_df = assign_lineages_de_novo(
        X_all,
        ensg_ids,
        obs_all,
        marker_rows,
        cfg["single_cell_processing"],
        random_seed=int(cfg["runtime_policy"]["random_seed"]),
    )

    out_disc = ROOT / "reports/discovery"
    out_disc.mkdir(parents=True, exist_ok=True)
    write_tsv(out_disc / "GSE179640_library_qc.tsv", lib_qc_rows, list(lib_qc_rows[0].keys()))
    cluster_df.to_csv(out_disc / "GSE179640_cluster_annotations.tsv", sep="\t", index=False)

    # Cell counts by donor x site x lineage
    dsl_counts = (
        annot_obs.groupby(["donor_id", "site_code", "lineage"], sort=True)
        .agg(n_cells=("cell_id", "count"), total_umis=("total_counts", "sum"))
        .reset_index()
    )
    dsl_counts.to_csv(out_disc / "GSE179640_donor_site_lineage_counts.tsv", sep="\t", index=False)

    # 4. Exact donor x site x lineage pseudobulk aggregation + conservation check
    valid_lineages = list(cfg["single_cell_processing"]["de_novo_lineage_annotation"]["lineages"])
    min_cells_pb = int(cfg["single_cell_processing"]["pseudobulk"]["min_cells_per_sample"])
    pb_by_lin, conservation_audit = aggregate_pseudobulk(
        X_all, ensg_ids, annot_obs, valid_lineages=valid_lineages, min_cells=min_cells_pb
    )

    # Save pseudobulk matrices in data/interim & sample ledger in reports/discovery
    interim_pb = ROOT / "data/interim/GSE179640_pseudobulk"
    interim_pb.mkdir(parents=True, exist_ok=True)
    pb_meta_all = []
    for lin, (mat, mdf) in pb_by_lin.items():
        np.savez_compressed(interim_pb / f"{lin}_counts.npz", counts=mat, ensg_ids=np.array(ensg_ids))
        pb_meta_all.append(mdf)
    pb_meta_df = pd.concat(pb_meta_all, ignore_index=True)
    pb_meta_df.to_csv(out_disc / "GSE179640_pseudobulk_samples.tsv", sep="\t", index=False)

    # Canonical marker purity audit across pseudobulk lineages (donor-equal-weighted mean CPM)
    gene_idx = {g: i for i, g in enumerate(ensg_ids)}
    purity_rows = []
    for m in marker_rows:
        gi = gene_idx.get(m["ensembl_gene_id"])
        if gi is None:
            continue
        row = {"marker_lineage": m["lineage"], "symbol": m["symbol"], "ensembl_gene_id": m["ensembl_gene_id"]}
        for lin in valid_lineages:
            if lin not in pb_by_lin:
                row[f"mean_cpm_{lin}"] = 0.0
                continue
            mat, _ = pb_by_lin[lin]
            cpm = mat[gi, :] / np.maximum(mat.sum(axis=0), 1.0) * 1e6
            row[f"mean_cpm_{lin}"] = round(float(cpm.mean()), 2)
        purity_rows.append(row)
    write_tsv(out_disc / "GSE179640_marker_purity_audit.tsv", purity_rows, list(purity_rows[0].keys()))

    # 5. Primary Stromal Fibroblast Paired Contrasts: (1) EcP vs EuE (n=8), (2) EcP vs AdjP (n=6)
    stromal_counts, stromal_meta = pb_by_lin["stromal_fibroblast"]
    paired_eue_ecp_sub = stromal_meta[stromal_meta["site_code"].isin(["EuE", "EcP"])]
    by_d_sites = paired_eue_ecp_sub.groupby("donor_id")["site_code"].apply(set)
    paired_8_donors = sorted(d for d, s in by_d_sites.items() if {"EuE", "EcP"} <= s)
    paired_8_cols = stromal_meta[
        stromal_meta["donor_id"].isin(paired_8_donors) & stromal_meta["site_code"].isin(["EuE", "EcP"])
    ].index.values

    # Expressed gene filter on paired EuE/EcP stromal pseudobulks: CPM >= 1.0 in >= 25% of samples
    C_p8 = stromal_counts[:, paired_8_cols]
    cpm_p8 = C_p8 / np.maximum(C_p8.sum(axis=0, keepdims=True), 1.0) * 1e6
    min_cpm = float(cfg["ranking_methods"]["gene_universe_rule"]["min_cpm"])
    min_frac = float(cfg["ranking_methods"]["gene_universe_rule"]["min_sample_fraction_in_paired_EuE_EcP"])
    expr_mask = (cpm_p8 >= min_cpm).mean(axis=1) >= min_frac
    expr_genes = [g for g, k in zip(ensg_ids, expr_mask) if k]
    expr_pos = {g: i for i, g in enumerate(expr_genes)}

    fit_eue = fit_paired_donor_voom(stromal_counts, stromal_meta, site_case="EcP", site_ref="EuE", keep_mask=expr_mask)
    fit_adjp = fit_paired_donor_voom(stromal_counts, stromal_meta, site_case="EcP", site_ref="AdjP", keep_mask=expr_mask)

    # Mesothelial ambient host RNA audit: mean log2(CPM + 1) in mesothelial pseudobulk vs EcP stromal pseudobulk
    ecp_stromal_cols = stromal_meta[stromal_meta["site_code"] == "EcP"].index.values
    ecp_stromal_cpm = stromal_counts[expr_mask][:, ecp_stromal_cols] / np.maximum(
        stromal_counts[:, ecp_stromal_cols].sum(axis=0, keepdims=True), 1.0
    ) * 1e6
    mean_log2cpm_ecp_stroma = np.log2(ecp_stromal_cpm + 1.0).mean(axis=1)

    if "mesothelial" in pb_by_lin:
        meso_counts, _ = pb_by_lin["mesothelial"]
        meso_cpm = meso_counts[expr_mask, :] / np.maximum(meso_counts.sum(axis=0, keepdims=True), 1.0) * 1e6
        mean_log2cpm_meso = np.log2(meso_cpm + 1.0).mean(axis=1)
    else:
        mean_log2cpm_meso = np.zeros(len(expr_genes), dtype=np.float64)

    delta_meso = mean_log2cpm_meso - mean_log2cpm_ecp_stroma
    w_ambient = 1.0 / (1.0 + np.exp(np.clip(delta_meso - 2.0, -30.0, 30.0)))

    # 6. Niche Partner Activity across Endothelial & Myeloid/Macrophage Lesion Compartments
    niche_compartments = cfg["ranking_methods"]["niche_partner_compartments"]
    min_partner_log2cpm1 = float(cfg["ranking_methods"]["methods"]["M2"]["niche_partner_min_log2cpm1"])
    comp_stats = {}
    for comp in niche_compartments:
        if comp not in pb_by_lin:
            continue
        c_mat, c_meta = pb_by_lin[comp]
        ecp_cols = c_meta[c_meta["site_code"] == "EcP"].index.values
        if len(ecp_cols) == 0:
            continue
        c_cpm_ecp = c_mat[:, ecp_cols] / np.maximum(c_mat[:, ecp_cols].sum(axis=0, keepdims=True), 1.0) * 1e6
        mean_l2_ecp = np.log2(c_cpm_ecp + 1.0).mean(axis=1)

        # Fit paired EcP vs EuE if >= 4 pairs, else donor-blocked GLS on EuE vs EcP
        sub_m = c_meta[c_meta["site_code"].isin(["EuE", "EcP"])].copy()
        by_d = sub_m.groupby("donor_id")["site_code"].apply(set)
        p_donors = sorted(d for d, s in by_d.items() if {"EuE", "EcP"} <= s)
        c_expr_mask = (c_mat[:, sub_m.index.values] / np.maximum(c_mat[:, sub_m.index.values].sum(axis=0, keepdims=True), 1.0) * 1e6 >= 1.0).mean(axis=1) >= 0.20
        if len(p_donors) >= 4:
            c_fit = fit_paired_donor_voom(c_mat, c_meta, site_case="EcP", site_ref="EuE", keep_mask=c_expr_mask)
        else:
            C_sub = c_mat[c_expr_mask][:, sub_m.index.values]
            l2_sub = np.log2(C_sub / np.maximum(c_mat[:, sub_m.index.values].sum(axis=0, keepdims=True), 1.0) * 1e6 + 1.0)
            X_sub = np.column_stack([np.ones(len(sub_m)), (sub_m["site_code"].values == "EcP").astype(float)])
            c_fit = limma_duplicate_correlation_gls(l2_sub, X_sub, sub_m["donor_id"].tolist(), coef_idx=1)

        c_genes = [g for g, k in zip(ensg_ids, c_expr_mask) if k]
        z_map = {g: float(z) for g, z in zip(c_genes, c_fit["z"])}
        q_map = {g: float(q) for g, q in zip(c_genes, c_fit["q"])}
        lfc_map = {g: float(b) for g, b in zip(c_genes, c_fit["logFC"])}
        comp_stats[comp] = {
            "mean_log2cpm1_ecp": {g: float(v) for g, v in zip(ensg_ids, mean_l2_ecp)},
            "z_ecp_vs_eue": z_map,
            "q_ecp_vs_eue": q_map,
            "logfc_ecp_vs_eue": lfc_map,
            "n_ecp_donors": int(len(ecp_cols)),
            "sig_q05_genes": int(sum(1 for q in c_fit["q"] if q < 0.05)),
        }

    # Also save compartment-level DE tables for Stage 3 genetic enrichment across compartments!
    for comp_name in ["stromal_fibroblast", *niche_compartments, "epithelial", "perivascular_mural"]:
        if comp_name not in pb_by_lin:
            continue
        c_mat, c_meta = pb_by_lin[comp_name]
        sub_m = c_meta[c_meta["site_code"].isin(["EuE", "EcP"])]
        by_d = sub_m.groupby("donor_id")["site_code"].apply(set)
        if sum(1 for s in by_d.tolist() if {"EuE", "EcP"} <= s) >= 4:
            cmask = (c_mat[:, sub_m.index.values] / np.maximum(c_mat[:, sub_m.index.values].sum(axis=0, keepdims=True), 1.0) * 1e6 >= 1.0).mean(axis=1) >= 0.25
            cfit = fit_paired_donor_voom(c_mat, c_meta, site_case="EcP", site_ref="EuE", keep_mask=cmask)
            cgenes = [g for g, k in zip(ensg_ids, cmask) if k]
            crows = [
                {
                    "ensembl_gene_id": g,
                    "symbol": hgnc["ensg_to_symbol"].get(g, ""),
                    "logFC_EcP_vs_EuE": round(float(cfit["logFC"][i]), 6),
                    "t": round(float(cfit["t"][i]), 6),
                    "z": round(float(cfit["z"][i]), 6),
                    "p": float(cfit["p"][i]),
                    "q_bh": float(cfit["q"][i]),
                }
                for i, g in enumerate(cgenes)
            ]
            write_tsv(out_disc / f"GSE179640_DE_{comp_name}_EcP_vs_EuE.tsv", crows, list(crows[0].keys()))

    # Evaluate curated L-R niche partner scores for every stromal gene
    niche_lr_by_gene = defaultdict(lambda: {"score": 0.0, "best_pair": "", "best_compartment": ""})
    lr_audit_rows = []
    stromal_ecp_l2_map = {g: float(v) for g, v in zip(expr_genes, mean_log2cpm_ecp_stroma)}

    for lr in lr_rows:
        lig_id, rec_id = lr["ligand_ensembl_id"], lr["receptor_ensembl_id"]
        for stromal_role, s_gid, p_gid, s_sym, p_sym in (
            ("ligand", lig_id, rec_id, lr["ligand_symbol"], lr["receptor_symbol"]),
            ("receptor", rec_id, lig_id, lr["receptor_symbol"], lr["ligand_symbol"]),
        ):
            if s_gid not in expr_pos:
                continue
            for comp, cst in comp_stats.items():
                p_l2 = cst["mean_log2cpm1_ecp"].get(p_gid, 0.0)
                if p_l2 < min_partner_log2cpm1:
                    continue
                p_z = cst["z_ecp_vs_eue"].get(p_gid, 0.0)
                act = min(1.0, max(0.0, p_l2 / 8.0)) * (1.0 / (1.0 + math.exp(-0.5 * max(-10.0, min(10.0, p_z)))))
                lr_audit_rows.append({
                    "pair_id": lr["pair_id"],
                    "pathway_group": lr["pathway_group"],
                    "stromal_gene_id": s_gid,
                    "stromal_symbol": s_sym,
                    "stromal_role": stromal_role,
                    "partner_gene_id": p_gid,
                    "partner_symbol": p_sym,
                    "niche_compartment": comp,
                    "stromal_ecp_log2cpm1": round(stromal_ecp_l2_map.get(s_gid, 0.0), 4),
                    "partner_ecp_log2cpm1": round(p_l2, 4),
                    "partner_z_ecp_vs_eue": round(p_z, 4),
                    "niche_activity_score": round(act, 6),
                })
                if act > niche_lr_by_gene[s_gid]["score"]:
                    niche_lr_by_gene[s_gid] = {
                        "score": act,
                        "best_pair": lr["pair_id"],
                        "best_compartment": comp,
                    }

    lr_audit_rows.sort(key=lambda r: -r["niche_activity_score"])
    if lr_audit_rows:
        write_tsv(out_disc / "GSE179640_niche_lr_scores.tsv", lr_audit_rows, list(lr_audit_rows[0].keys()))

    # 7. Shared Comparative Universe U on STRING v12.5 direct functional association network
    net_path = ROOT / cfg["ranking_methods"]["gene_universe_rule"]["network_reference"]
    expr_set = set(expr_genes)
    induced_edges = []
    with gzip.open(net_path, "rt", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            ga, gb, w = r["gene1"], r["gene2"], float(r["weight"])
            if ga in expr_set and gb in expr_set:
                induced_edges.append((ga, gb, w))

    universe = sorted({g for ga, gb, _ in induced_edges for g in (ga, gb)})
    u_idx = {g: i for i, g in enumerate(universe)}
    left_idx = np.array([u_idx[ga] for ga, gb, _ in induced_edges], dtype=np.int64)
    right_idx = np.array([u_idx[gb] for ga, gb, _ in induced_edges], dtype=np.int64)
    edge_w = np.array([w for _, _, w in induced_edges], dtype=np.float64)
    w_deg = np.bincount(left_idx, weights=edge_w, minlength=len(universe)) + np.bincount(
        right_idx, weights=edge_w, minlength=len(universe)
    )

    # Extract Discovery vectors aligned to shared universe U
    u_pos_in_expr = np.array([expr_pos[g] for g in universe], dtype=np.int64)
    lfc_eue = fit_eue["logFC"][u_pos_in_expr]
    t_eue = fit_eue["t"][u_pos_in_expr]
    z_eue = fit_eue["z"][u_pos_in_expr]
    p_eue = fit_eue["p"][u_pos_in_expr]
    q_eue = fit_eue["q"][u_pos_in_expr]

    lfc_adjp = fit_adjp["logFC"][u_pos_in_expr]
    t_adjp = fit_adjp["t"][u_pos_in_expr]
    z_adjp = fit_adjp["z"][u_pos_in_expr]
    p_adjp = fit_adjp["p"][u_pos_in_expr]
    q_adjp = fit_adjp["q"][u_pos_in_expr]

    w_amb_u = w_ambient[u_pos_in_expr]
    delta_meso_u = delta_meso[u_pos_in_expr]

    disc_sign = np.where(lfc_eue >= 0.0, 1, -1)
    z1_abs = np.abs(z_eue)
    z2_aligned = np.sign(t_eue) * z_adjp

    # Method 1: B0 (Naive Paired Lesion DE: EcP vs EuE)
    score_B0 = z1_abs.copy()

    # Method 2: G1 (v1 Historical Baseline: Global STRING RWR on abs(t_eue))
    score_G1, rwr_diag = rwr_propagate(np.abs(t_eue), left_idx, right_idx, edge_w, restart_prob=0.5, tol=1e-10)

    # Method 3: M1 (Host-Decontaminated Dual-Contrast Lesion Score)
    m1_core = w_amb_u * np.minimum(z1_abs, np.maximum(0.0, z2_aligned))
    tie_term = 1e-4 * (z1_abs + z2_aligned)
    score_M1 = m1_core + tie_term

    # Method 4: M2 (Genetics & Niche Crosstalk-Anchored Target Score)
    alpha_gen = float(cfg["ranking_methods"]["methods"]["M2"]["alpha_genetics"])
    alpha_lr = float(cfg["ranking_methods"]["methods"]["M2"]["alpha_niche_lr"])
    g_prior_u = np.array(
        [float(genetics_by_ensg[g]["genetic_prior_score"]) if g in genetics_by_ensg else 0.0 for g in universe],
        dtype=np.float64,
    )
    n_lr_u = np.array([float(niche_lr_by_gene[g]["score"]) for g in universe], dtype=np.float64)
    score_M2 = m1_core * (1.0 + alpha_gen * g_prior_u + alpha_lr * n_lr_u) + tie_term

    scores_by_method = {"B0": score_B0, "G1": score_G1, "M1": score_M1, "M2": score_M2}
    ranks_by_method = {}
    for mname, sc_vec in scores_by_method.items():
        order = sorted(range(len(universe)), key=lambda i: (-float(sc_vec[i]), universe[i]))
        r_arr = np.empty(len(universe), dtype=np.int64)
        for rk, idx_u in enumerate(order, 1):
            r_arr[idx_u] = rk
        ranks_by_method[mname] = r_arr

    # 8. Write full universe ranking table + frozen shortlist manifest (K=20 and K=50)
    full_rows = []
    shortlist_rows = []
    for i, gid in enumerate(universe):
        sym = hgnc["ensg_to_symbol"].get(gid, "")
        gen_rec = genetics_by_ensg.get(gid, {})
        row = {
            "ensembl_gene_id": gid,
            "symbol": sym,
            "discovery_sign": int(disc_sign[i]),
            "logFC_EcP_vs_EuE": round(float(lfc_eue[i]), 6),
            "t_EcP_vs_EuE": round(float(t_eue[i]), 6),
            "z_EcP_vs_EuE": round(float(z_eue[i]), 6),
            "q_EcP_vs_EuE": float(q_eue[i]),
            "logFC_EcP_vs_AdjP": round(float(lfc_adjp[i]), 6),
            "t_EcP_vs_AdjP": round(float(t_adjp[i]), 6),
            "z_EcP_vs_AdjP": round(float(z_adjp[i]), 6),
            "q_EcP_vs_AdjP": float(q_adjp[i]),
            "delta_meso_minus_stroma_log2cpm": round(float(delta_meso_u[i]), 4),
            "w_ambient": round(float(w_amb_u[i]), 6),
            "m1_core": round(float(m1_core[i]), 6),
            "genetic_prior_score": round(float(g_prior_u[i]), 6),
            "niche_lr_score": round(float(n_lr_u[i]), 6),
            "best_lr_pair": niche_lr_by_gene[gid]["best_pair"],
            "best_lr_compartment": niche_lr_by_gene[gid]["best_compartment"],
            "string_weighted_degree": round(float(w_deg[i]), 4),
            "score_B0": float(score_B0[i]),
            "rank_B0": int(ranks_by_method["B0"][i]),
            "score_G1": float(score_G1[i]),
            "rank_G1": int(ranks_by_method["G1"][i]),
            "score_M1": float(score_M1[i]),
            "rank_M1": int(ranks_by_method["M1"][i]),
            "score_M2": float(score_M2[i]),
            "rank_M2": int(ranks_by_method["M2"][i]),
            "in_B0_top20": str(ranks_by_method["B0"][i] <= 20),
            "in_B0_top50": str(ranks_by_method["B0"][i] <= 50),
            "in_G1_top20": str(ranks_by_method["G1"][i] <= 20),
            "in_G1_top50": str(ranks_by_method["G1"][i] <= 50),
            "in_M1_top20": str(ranks_by_method["M1"][i] <= 20),
            "in_M1_top50": str(ranks_by_method["M1"][i] <= 50),
            "in_M2_top20": str(ranks_by_method["M2"][i] <= 20),
            "in_M2_top50": str(ranks_by_method["M2"][i] <= 50),
            "has_human_genetics": gen_rec.get("has_human_genetic_support", "False"),
            "tractable_sm_or_ab": str(
                gen_rec.get("tractable_sm_any") == "True" or gen_rec.get("tractable_ab_any") == "True"
            ),
        }
        full_rows.append(row)
        if min(row["rank_B0"], row["rank_G1"], row["rank_M1"], row["rank_M2"]) <= 50:
            shortlist_rows.append(row)

    shortlist_rows.sort(key=lambda r: (min(r["rank_M2"], r["rank_M1"], r["rank_B0"], r["rank_G1"]), r["ensembl_gene_id"]))

    full_gz_path = out_disc / "GSE179640_stromal_rankings_full.tsv.gz"
    with full_gz_path.open("wb") as raw_f, gzip.GzipFile(filename="", mode="wb", fileobj=raw_f, mtime=0) as gz:
        buf = pd.DataFrame(full_rows).to_csv(sep="\t", index=False, lineterminator="\n").encode("utf-8")
        gz.write(buf)

    frozen_cand_path = ROOT / "manifests/frozen_candidates_v2.tsv"
    write_tsv(frozen_cand_path, shortlist_rows, list(shortlist_rows[0].keys()))

    # Summarize discovery statistics and inspect host-confounding diagnostics in B0 vs M1/M2
    def summarize_method_top(mname: str, k: int):
        sub = [r for r in shortlist_rows if r[f"rank_{mname}"] <= k]
        sub.sort(key=lambda r: r[f"rank_{mname}"])
        return {
            "k": k,
            "top10_symbols": [r["symbol"] or r["ensembl_gene_id"] for r in sub[:10]],
            "dual_concordant_count": sum(1 for r in sub if r["m1_core"] > 0),
            "ambient_meso_enriched_gt2_count": sum(1 for r in sub if r["delta_meso_minus_stroma_log2cpm"] > 2.0),
            "human_genetics_supported_count": sum(1 for r in sub if r["has_human_genetics"] == "True"),
            "active_niche_lr_count": sum(1 for r in sub if r["niche_lr_score"] > 0),
            "tractable_sm_or_ab_count": sum(1 for r in sub if r["tractable_sm_or_ab"] == "True"),
            "median_string_degree": round(float(np.median([r["string_weighted_degree"] for r in sub])), 2),
        }

    disc_summary = {
        "schema_version": 2,
        "cohort": "GSE179640",
        "libraries_processed": len(lib_qc_rows),
        "raw_cells_total": int(sum(r["raw_cells"] for r in lib_qc_rows)),
        "qc_pass_cells_total": int(X_all.shape[0]),
        "lineage_cell_counts": {k: int(v) for k, v in Counter(annot_obs["lineage"]).items()},
        "pseudobulk_conservation": conservation_audit,
        "paired_stromal_models": {
            "EcP_vs_EuE": {
                "n_pairs": fit_eue["n_pairs"],
                "paired_donors": fit_eue["paired_donors"],
                "df_residual": fit_eue["df_residual"],
                "prior_df": round(fit_eue["prior_df"], 4),
                "sig_BH_q_lt_0.05_expressed": int((fit_eue["q"] < 0.05).sum()),
                "sig_BH_q_lt_0.05_universe": int((q_eue < 0.05).sum()),
            },
            "EcP_vs_AdjP": {
                "n_pairs": fit_adjp["n_pairs"],
                "paired_donors": fit_adjp["paired_donors"],
                "df_residual": fit_adjp["df_residual"],
                "prior_df": round(fit_adjp["prior_df"], 4),
                "sig_BH_q_lt_0.05_expressed": int((fit_adjp["q"] < 0.05).sum()),
                "sig_BH_q_lt_0.05_universe": int((q_adjp < 0.05).sum()),
            },
        },
        "gene_universe": {
            "expressed_stromal_genes": len(expr_genes),
            "shared_string_universe_genes": len(universe),
            "induced_string_edges": len(induced_edges),
            "dual_contrast_q05_both_genes": int(((q_eue < 0.05) & (q_adjp < 0.05) & (z2_aligned > 0)).sum()),
        },
        "rwr_propagation": rwr_diag,
        "method_shortlist_diagnostics": {
            m: {"top20": summarize_method_top(m, 20), "top50": summarize_method_top(m, 50)}
            for m in ("B0", "G1", "M1", "M2")
        },
    }
    disc_sum_path = out_disc / "discovery_summary_v2.json"
    disc_sum_path.write_text(json.dumps(disc_summary, indent=2) + "\n", encoding="utf-8")

    # 9. Seal discovery predictions BEFORE external data download/evaluation
    sealed_files = [
        "configs/protocol_v2.yaml",
        "receipts/protocol_freeze.json",
        "receipts/discovery_GSE179640_data.json",
        "manifests/frozen_candidates_v2.tsv",
        "reports/discovery/GSE179640_stromal_rankings_full.tsv.gz",
        "reports/discovery/GSE179640_niche_lr_scores.tsv",
        "reports/discovery/GSE179640_library_qc.tsv",
        "reports/discovery/GSE179640_cluster_annotations.tsv",
        "reports/discovery/GSE179640_donor_site_lineage_counts.tsv",
        "reports/discovery/GSE179640_pseudobulk_samples.tsv",
        "reports/discovery/GSE179640_marker_purity_audit.tsv",
        "reports/discovery/discovery_summary_v2.json",
        "scripts/endotarget_lib.py",
        "scripts/run_discovery_v2.py",
    ]
    seal = {
        "schema_version": 2,
        "phase": "STAGE2_PREDICTION_SEAL",
        "sealed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "protocol_freeze_sha256": freeze_sha,
        "discovery_data_receipt_sha256": sha256_file(data_rec_path),
        "attestations": {
            "external_GSE213216_downloaded_or_read": False,
            "external_GSE141549_expression_downloaded_or_read": False,
            "heca_annotations_used": False,
            "pseudobulk_umi_conservation_exact": True,
        },
        "sealed_artifacts": {
            rel: {"sha256": sha256_file(ROOT / rel), "bytes": (ROOT / rel).stat().st_size}
            for rel in sealed_files
        },
    }
    seal_path = ROOT / "receipts/prediction_seal.json"
    seal_path.write_text(json.dumps(seal, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "SEALED",
        "seal_sha256": sha256_file(seal_path),
        "qc_pass_cells": disc_summary["qc_pass_cells_total"],
        "universe_genes": len(universe),
        "paired_models": disc_summary["paired_stromal_models"],
        "shortlist_top20": {m: disc_summary["method_shortlist_diagnostics"][m]["top20"] for m in ("B0", "G1", "M1", "M2")},
    }, indent=2))


if __name__ == "__main__":
    main()
