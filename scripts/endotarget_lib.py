#!/usr/bin/env python3
"""Core statistical, single-cell de novo annotation, pseudobulk, limma-voom,
duplicate-correlation GLS, and STRING RWR library for EndoTarget v2.

Strictly enforces thread cap (<= 8 threads) and zero HECA annotation usage.
"""
from __future__ import annotations

import csv
import gzip
import math
import os
import pathlib
import re
import tarfile
from collections import defaultdict

os.environ["OMP_NUM_THREADS"] = "8"
os.environ["OPENBLAS_NUM_THREADS"] = "8"
os.environ["MKL_NUM_THREADS"] = "8"
os.environ["NUMEXPR_NUM_THREADS"] = "8"

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
from scipy import special, stats

sc.settings.n_jobs = 8
ENSG_RE = re.compile(r"^ENSG[0-9]{11}$")


# ---------------------------------------------------------------------------
# 1. Raw 10x Matrix Loaders (HDF5 and tar.gz MTX) -> canonical ENSG CSR counts
# ---------------------------------------------------------------------------
def strip_ensg_version(gid: str) -> str:
    g = str(gid).strip()
    if g.startswith("ENSG") and "." in g:
        g = g.split(".", 1)[0]
    return g


def load_10x_h5_canonical(path: pathlib.Path, sym_to_ensg: dict[str, str] | None = None):
    """Load 10x Cell Ranger filtered_feature_bc_matrix.h5 into (X_csr_int64, barcodes, ensg_ids, is_mt_raw_cell_pct)."""
    with h5py.File(path, "r") as f:
        grp = f["matrix"]
        data = grp["data"][:]
        indices = grp["indices"][:]
        indptr = grp["indptr"][:]
        shape = tuple(grp["shape"][:])  # (n_features, n_cells) in 10x h5
        barcodes = [b.decode("utf-8") if isinstance(b, bytes) else str(b) for b in grp["barcodes"][:]]
        feat = grp["features"]
        raw_ids = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in feat["id"][:]]
        raw_names = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in feat["name"][:]]
        if "feature_type" in feat:
            ftypes = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in feat["feature_type"][:]]
        else:
            ftypes = ["Gene Expression"] * len(raw_ids)

    # CSC over (features x cells) == CSR over (cells x features)!
    X_cells_by_feat = sp.csr_matrix((data.astype(np.int64), indices, indptr), shape=(shape[1], shape[0]))
    return _canonicalize_raw_counts(X_cells_by_feat, barcodes, raw_ids, raw_names, ftypes, sym_to_ensg)


def load_10x_targz_or_mtx_canonical(
    paths: list[pathlib.Path],
    sym_to_ensg: dict[str, str] | None = None,
):
    """Load 10x Cell Ranger output from either a single .tar.gz archive or separate barcodes/features/matrix.gz files."""
    if len(paths) == 1 and paths[0].name.endswith(".tar.gz"):
        with tarfile.open(paths[0], "r:gz") as tf:
            members = {m.name: m for m in tf.getmembers() if m.isfile()}
            # Prefer filtered_feature_bc_matrix if both raw and filtered exist inside tar.gz
            cand_mtx = [n for n in members if n.endswith("matrix.mtx.gz") or n.endswith("matrix.mtx")]
            filt_mtx = [n for n in cand_mtx if "raw_" not in n]
            mtx_name = sorted(filt_mtx or cand_mtx)[0]
            prefix = mtx_name[: mtx_name.rfind("matrix.mtx")]
            bc_candidates = [n for n in members if n.startswith(prefix) and "barcodes.tsv" in n]
            ft_candidates = [n for n in members if n.startswith(prefix) and ("features.tsv" in n or "genes.tsv" in n)]
            if not bc_candidates:
                bc_candidates = [n for n in members if "barcodes.tsv" in n and "raw_" not in n]
            if not ft_candidates:
                ft_candidates = [n for n in members if ("features.tsv" in n or "genes.tsv" in n) and "raw_" not in n]
            bc_name = sorted(bc_candidates)[0]
            ft_name = sorted(ft_candidates)[0]

            def read_lines_member(mname: str) -> list[str]:
                raw = tf.extractfile(members[mname]).read()
                if mname.endswith(".gz"):
                    raw = gzip.decompress(raw)
                return raw.decode("utf-8", errors="replace").splitlines()

            barcodes = [l.strip() for l in read_lines_member(bc_name) if l.strip()]
            ft_lines = [l.rstrip("\r\n").split("\t") for l in read_lines_member(ft_name) if l.strip()]
            raw_mtx = tf.extractfile(members[mtx_name]).read()
            if mtx_name.endswith(".gz"):
                raw_mtx = gzip.decompress(raw_mtx)
            import io
            mat_coo = scipy.io.mmread(io.BytesIO(raw_mtx))
    else:
        bc_p = next(p for p in paths if "barcodes.tsv" in p.name)
        ft_p = next(p for p in paths if "features.tsv" in p.name or "genes.tsv" in p.name)
        mt_p = next(p for p in paths if "matrix.mtx" in p.name)
        with gzip.open(bc_p, "rt", encoding="utf-8") as f:
            barcodes = [l.strip() for l in f if l.strip()]
        with gzip.open(ft_p, "rt", encoding="utf-8") as f:
            ft_lines = [l.rstrip("\r\n").split("\t") for l in f if l.strip()]
        with gzip.open(mt_p, "rb") as f:
            mat_coo = scipy.io.mmread(f)

    raw_ids = [r[0] for r in ft_lines]
    raw_names = [r[1] if len(r) > 1 else r[0] for r in ft_lines]
    ftypes = [r[2] if len(r) > 2 else "Gene Expression" for r in ft_lines]
    X_cells_by_feat = mat_coo.T.tocsr().astype(np.int64)
    return _canonicalize_raw_counts(X_cells_by_feat, barcodes, raw_ids, raw_names, ftypes, sym_to_ensg)


def _canonicalize_raw_counts(
    X_cells_by_feat: sp.csr_matrix,
    barcodes: list[str],
    raw_ids: list[str],
    raw_names: list[str],
    ftypes: list[str],
    sym_to_ensg: dict[str, str] | None,
):
    gex_mask = np.array([ft in ("Gene Expression", "") for ft in ftypes], dtype=bool)
    X_gex = X_cells_by_feat[:, gex_mask]
    names_gex = [n for n, m in zip(raw_names, gex_mask) if m]
    ids_gex = [i for i, m in zip(raw_ids, gex_mask) if m]

    # Compute QC statistics BEFORE dropping unmapped non-ENSG transcripts so MT% and total_counts are exact!
    total_counts = np.asarray(X_gex.sum(axis=1)).ravel().astype(np.int64)
    n_genes_by_counts = np.diff(X_gex.indptr).astype(np.int64)
    mt_mask = np.array([str(n).upper().startswith("MT-") for n in names_gex], dtype=bool)
    mt_counts = np.asarray(X_gex[:, mt_mask].sum(axis=1)).ravel().astype(np.int64) if mt_mask.any() else np.zeros(X_gex.shape[0], dtype=np.int64)
    pct_mt = np.where(total_counts > 0, 100.0 * mt_counts / np.maximum(total_counts, 1), 0.0)

    # Map features to canonical unversioned ENSG IDs and sum duplicate ENSG features if any
    ensg_list = []
    keep_cols = []
    for col_idx, (gid, gname) in enumerate(zip(ids_gex, names_gex)):
        cand = strip_ensg_version(gid)
        if not ENSG_RE.fullmatch(cand) and sym_to_ensg is not None:
            cand = sym_to_ensg.get(gname.strip(), "")
        if ENSG_RE.fullmatch(cand):
            keep_cols.append(col_idx)
            ensg_list.append(cand)

    X_sub = X_gex[:, keep_cols]
    unique_ensg = sorted(set(ensg_list))
    if len(unique_ensg) == len(ensg_list) and ensg_list == unique_ensg:
        X_canon = X_sub
    else:
        pos_map = {g: i for i, g in enumerate(unique_ensg)}
        col_targets = np.array([pos_map[g] for g in ensg_list], dtype=np.int64)
        # Sparse projection matrix (n_sub x n_unique) to collapse duplicates & sort deterministically
        proj = sp.csr_matrix(
            (np.ones(len(keep_cols), dtype=np.int64), (np.arange(len(keep_cols), dtype=np.int64), col_targets)),
            shape=(len(keep_cols), len(unique_ensg)),
        )
        X_canon = (X_sub @ proj).tocsr()

    qc_df = pd.DataFrame({
        "barcode": barcodes,
        "total_counts": total_counts,
        "n_genes_by_counts": n_genes_by_counts,
        "pct_counts_mt": pct_mt,
    })
    return X_canon, unique_ensg, qc_df


# ---------------------------------------------------------------------------
# 2. De Novo Canonical Lineage Assignment (Zero HECA Labels)
# ---------------------------------------------------------------------------
def assign_lineages_de_novo(
    X_counts: sp.csr_matrix,
    ensg_ids: list[str],
    obs_df: pd.DataFrame,
    marker_rows: list[dict],
    cfg_sc: dict,
    random_seed: int = 20261120,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Perform de novo Leiden clustering + canonical marker z-score lineage assignment + cell purity guard.
    Assumes X_counts already restricted to QC-passing cells.
    Returns (annotated_obs_df, cluster_summary_df).
    """
    n_cells, n_genes = X_counts.shape
    adata = ad.AnnData(X=X_counts.astype(np.float32), obs=obs_df.copy())
    adata.var_names = pd.Index(ensg_ids)

    sc.pp.normalize_total(adata, target_sum=float(cfg_sc["de_novo_lineage_annotation"]["normalization_target_sum"]))
    sc.pp.log1p(adata)

    # Extract dense log1p expression ONLY for the 56 canonical markers (~30-50 MB even on 250k cells!)
    gene_to_idx = {g: i for i, g in enumerate(ensg_ids)}
    lineages = list(cfg_sc["de_novo_lineage_annotation"]["lineages"])
    markers_by_lin = defaultdict(list)
    for r in marker_rows:
        if r["ensembl_gene_id"] in gene_to_idx:
            markers_by_lin[r["lineage"]].append(r["ensembl_gene_id"])

    present_markers = sorted({g for gs in markers_by_lin.values() for g in gs})
    m_indices = [gene_to_idx[g] for g in present_markers]
    M_dense = adata.X[:, m_indices].toarray() if sp.issparse(adata.X) else np.asarray(adata.X[:, m_indices])
    m_pos = {g: i for i, g in enumerate(present_markers)}

    # Standardize each canonical marker across cells, clipped to [-3, +6]
    lo_clip, hi_clip = cfg_sc["de_novo_lineage_annotation"]["gene_zscore_clip"]
    means = M_dense.mean(axis=0, keepdims=True)
    stds = M_dense.std(axis=0, keepdims=True)
    stds = np.where(stds < 1e-6, 1.0, stds)
    Z_markers = np.clip((M_dense - means) / stds, lo_clip, hi_clip)

    cell_lin_scores = np.zeros((n_cells, len(lineages)), dtype=np.float32)
    for j, lin in enumerate(lineages):
        cols = [m_pos[g] for g in markers_by_lin[lin]]
        if cols:
            cell_lin_scores[:, j] = Z_markers[:, cols].mean(axis=1)

    # Unsupervised HVG -> PCA -> Neighbors -> Leiden (igraph backend, deterministic seed)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=min(int(cfg_sc["de_novo_lineage_annotation"]["hvg_n_top_genes"]), n_genes),
        flavor="seurat",
    )
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10.0)
    n_comps = min(int(cfg_sc["de_novo_lineage_annotation"]["pca_n_comps"]), adata_hvg.shape[1] - 1, adata_hvg.shape[0] - 1)
    sc.tl.pca(adata_hvg, n_comps=n_comps, random_state=random_seed, svd_solver="arpack")
    sc.pp.neighbors(
        adata_hvg,
        n_neighbors=int(cfg_sc["de_novo_lineage_annotation"]["neighbors_n"]),
        n_pcs=n_comps,
        random_state=random_seed,
    )
    sc.tl.leiden(
        adata_hvg,
        resolution=float(cfg_sc["de_novo_lineage_annotation"]["leiden_resolution"]),
        random_state=random_seed,
        flavor="igraph",
        n_iterations=2,
    )
    clusters = adata_hvg.obs["leiden"].astype(str).values

    min_cl_z = float(cfg_sc["de_novo_lineage_annotation"]["min_cluster_top_z"])
    max_margin = float(cfg_sc["de_novo_lineage_annotation"]["cell_cross_lineage_margin_max"])

    cluster_rows = []
    cluster_to_lin = {}
    cluster_to_idx = {}
    for cl in sorted(set(clusters), key=lambda x: int(x) if x.isdigit() else x):
        mask = clusters == cl
        cl_means = cell_lin_scores[mask].mean(axis=0)
        best_j = int(np.argmax(cl_means))
        best_z = float(cl_means[best_j])
        assigned = lineages[best_j] if best_z >= min_cl_z else "unassigned_low_score"
        cluster_to_lin[cl] = assigned
        cluster_to_idx[cl] = best_j
        row = {"cluster": cl, "n_cells": int(mask.sum()), "cluster_lineage": assigned, "top_z": round(best_z, 4)}
        for j, lin in enumerate(lineages):
            row[f"z_{lin}"] = round(float(cl_means[j]), 4)
        cluster_rows.append(row)

    # Cell-level purity guard: reject cells whose individual best lineage score beats the cluster-assigned lineage by > max_margin
    final_lineage = []
    cell_top_j = np.argmax(cell_lin_scores, axis=1)
    cell_top_z = cell_lin_scores[np.arange(n_cells), cell_top_j]
    for i in range(n_cells):
        cl = clusters[i]
        cl_lin = cluster_to_lin[cl]
        if cl_lin == "unassigned_low_score":
            final_lineage.append("unassigned_low_score")
            continue
        cl_j = cluster_to_idx[cl]
        own_cl_z = float(cell_lin_scores[i, cl_j])
        if cell_top_j[i] != cl_j and (float(cell_top_z[i]) - own_cl_z) > max_margin:
            final_lineage.append("ambiguous_doublet")
        else:
            final_lineage.append(cl_lin)

    out_obs = obs_df.copy()
    out_obs["leiden_cluster"] = clusters
    out_obs["cluster_lineage"] = [cluster_to_lin[c] for c in clusters]
    out_obs["lineage"] = final_lineage
    for j, lin in enumerate(lineages):
        out_obs[f"score_{lin}"] = cell_lin_scores[:, j]

    return out_obs, pd.DataFrame(cluster_rows)


# ---------------------------------------------------------------------------
# 3. Exact Donor x Site x Lineage Pseudobulk Aggregation
# ---------------------------------------------------------------------------
def aggregate_pseudobulk(
    X_counts: sp.csr_matrix,
    ensg_ids: list[str],
    obs_df: pd.DataFrame,
    valid_lineages: list[str],
    min_cells: int = 20,
) -> tuple[dict[str, tuple[np.ndarray, pd.DataFrame]], dict]:
    """Sum exact raw integer counts per (donor_id, site_code, lineage).
    Returns:
      by_lineage[lineage] = (counts_genes_x_samples_int64, sample_meta_df)
      conservation_audit = dict verifying exact integer count conservation!
    """
    total_raw_sum = int(X_counts.sum())
    kept_sum = 0
    dropped_low_cells_sum = 0
    excluded_lineage_sum = 0

    valid_set = set(valid_lineages)
    excl_mask = ~obs_df["lineage"].isin(valid_set).values
    if excl_mask.any():
        excluded_lineage_sum = int(X_counts[excl_mask].sum())

    by_lineage = {}
    for lin in valid_lineages:
        lin_mask = (obs_df["lineage"] == lin).values
        if not lin_mask.any():
            continue
        sub_obs = obs_df.loc[lin_mask]
        sub_X = X_counts[lin_mask]
        groups = sub_obs.groupby(["donor_id", "site_code"], sort=True).indices
        cols = []
        meta_rows = []
        for (donor_id, site_code), idx in sorted(groups.items()):
            c_sum = np.asarray(sub_X[idx].sum(axis=0)).ravel().astype(np.int64)
            if len(idx) < min_cells:
                dropped_low_cells_sum += int(c_sum.sum())
                continue
            kept_sum += int(c_sum.sum())
            cols.append(c_sum)
            meta_rows.append({
                "sample_id": f"{donor_id}__{site_code}__{lin}",
                "donor_id": donor_id,
                "site_code": site_code,
                "lineage": lin,
                "n_cells": int(len(idx)),
                "library_size": int(c_sum.sum()),
            })
        if cols:
            mat = np.column_stack(cols)  # genes x samples
            by_lineage[lin] = (mat, pd.DataFrame(meta_rows))

    assert total_raw_sum == kept_sum + dropped_low_cells_sum + excluded_lineage_sum, (
        f"Pseudobulk count conservation violated: {total_raw_sum} != {kept_sum} + {dropped_low_cells_sum} + {excluded_lineage_sum}"
    )
    audit = {
        "total_qc_passing_umi_sum": total_raw_sum,
        "retained_pseudobulk_umi_sum": kept_sum,
        "dropped_below_min_cells_umi_sum": dropped_low_cells_sum,
        "excluded_ambiguous_or_unassigned_umi_sum": excluded_lineage_sum,
        "exact_conservation": True,
    }
    return by_lineage, audit


# ---------------------------------------------------------------------------
# 4. TMM, limma-voom, Empirical Bayes Moderation, and Donor-Blocked GLS
# ---------------------------------------------------------------------------
def calc_tmm_norm_factors(counts: np.ndarray, logratio_trim: float = 0.30, sum_trim: float = 0.05) -> np.ndarray:
    """Robinson & Oshlack (2010) TMM normalization factors on genes x samples integer count matrix."""
    counts = np.asarray(counts, dtype=np.float64)
    lib_sizes = counts.sum(axis=0)
    if np.any(lib_sizes <= 0):
        raise ValueError("Zero library size in pseudobulk counts")
    # Reference column: sample whose 75th-percentile scaled count is closest to mean 75th percentile
    positive_genes = counts.sum(axis=1) > 0
    C = counts[positive_genes]
    q75 = np.percentile(C / lib_sizes[None, :], 75, axis=0)
    ref_idx = int(np.argmin(np.abs(q75 - q75.mean())))

    ref = C[:, ref_idx]
    n_ref = lib_sizes[ref_idx]
    factors = np.ones(C.shape[1], dtype=np.float64)

    for j in range(C.shape[1]):
        if j == ref_idx:
            continue
        obs = C[:, j]
        n_obs = lib_sizes[j]
        both_pos = (obs > 0) & (ref > 0)
        if both_pos.sum() < 20:
            factors[j] = 1.0
            continue
        o = obs[both_pos]
        r = ref[both_pos]
        logR = np.log2((o / n_obs) / (r / n_ref))
        absE = 0.5 * np.log2((o / n_obs) * (r / n_ref))
        v = (n_obs - o) / (n_obs * o) + (n_ref - r) / (n_ref * r)

        fin = np.isfinite(logR) & np.isfinite(absE) & (v > 0)
        logR, absE, v = logR[fin], absE[fin], v[fin]
        n_k = len(logR)
        if n_k < 10:
            factors[j] = 1.0
            continue
        lo_L = int(math.floor(n_k * logratio_trim)) + 1
        hi_L = n_k + 1 - lo_L
        lo_S = int(math.floor(n_k * sum_trim)) + 1
        hi_S = n_k + 1 - lo_S
        r_L = stats.rankdata(logR, method="ordinal")
        r_S = stats.rankdata(absE, method="ordinal")
        keep = (r_L >= lo_L) & (r_L <= hi_L) & (r_S >= lo_S) & (r_S <= hi_S)
        if keep.sum() < 5:
            factors[j] = 1.0
        else:
            w = 1.0 / v[keep]
            factors[j] = 2.0 ** float(np.sum(w * logR[keep]) / np.sum(w))

    # Geometric mean = 1 scaling (edgeR convention)
    factors = factors / np.exp(np.mean(np.log(factors)))
    return factors


def trigamma_inverse(y: float) -> float:
    """Exact Smyth (2004) / limma::trigammaInverse Newton solver."""
    if y > 1e7:
        return 1.0 / math.sqrt(y)
    if y < 1e-6:
        return 1.0 / y
    x = 0.5 + 1.0 / y
    for _ in range(50):
        tri = float(special.polygamma(1, x))
        dif = tri * (1.0 - tri / y) / float(special.polygamma(2, x))
        x = x + dif
        if abs(-dif / x) < 1e-8:
            break
    return float(x)


def squeeze_var(s2: np.ndarray, df: float) -> tuple[np.ndarray, float, float]:
    """Smyth (2004) empirical Bayes variance moderation."""
    s2 = np.asarray(s2, dtype=np.float64)
    ok = np.isfinite(s2) & (s2 > 1e-15)
    z = np.log(s2[ok])
    e = z - float(special.digamma(df / 2.0)) + math.log(df / 2.0)
    emean = float(e.mean())
    evar = max(float(np.var(e, ddof=1)) - float(special.polygamma(1, df / 2.0)), 1e-8)
    x = trigamma_inverse(evar)
    d0 = 2.0 * x
    s0 = math.exp(emean + float(special.digamma(d0 / 2.0)) - math.log(d0 / 2.0))
    s2_post = (d0 * s0 + df * s2) / (d0 + df)
    return s2_post, float(d0), float(s0)


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=np.float64)
    m = p.size
    if m == 0:
        return np.array([], dtype=np.float64)
    order = np.argsort(p)
    ranked = p[order] * m / np.arange(1, m + 1, dtype=np.float64)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m, dtype=np.float64)
    out[order] = np.clip(q, 0.0, 1.0)
    return out


def signed_z_from_t_df(t_vals: np.ndarray, df_total: float) -> np.ndarray:
    """Convert Student-t statistics with df_total to exact two-sided normal z-equivalents
    preserving sign and avoiding precision underflow in extreme tails."""
    t_vals = np.asarray(t_vals, dtype=np.float64)
    log_sf = stats.t.logsf(np.abs(t_vals), df=df_total)
    # stats.norm.isf(exp(log_sf)): use ndtri via log survival function when p is tiny
    # For moderate |t|, isf(sf) is exact; for extreme |t|, asymptotic log_sf matching via norm.isf
    sf = np.exp(log_sf)
    z_abs = np.where(
        sf > 1e-300,
        stats.norm.isf(np.clip(sf, 1e-300, 0.5)),
        np.sqrt(np.maximum(0.0, -2.0 * log_sf - np.log(2.0 * np.pi))),
    )
    return np.sign(t_vals) * z_abs


def _lowess_1d_smooth(x: np.ndarray, y: np.ndarray, n_grid: int = 200, span: float = 0.5) -> tuple[np.ndarray, np.ndarray]:
    """Deterministic tricube-weighted local linear smooth on an evenly spaced x-grid for voom mean-variance trend."""
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    grid = np.linspace(float(xs[0]), float(xs[-1]), n_grid)
    k = max(10, int(round(span * len(xs))))
    pred = np.empty(n_grid, dtype=np.float64)
    for idx_g, gx in enumerate(grid):
        dists = np.abs(xs - gx)
        # Partition top-k nearest points in O(N)
        part = np.argpartition(dists, k - 1)[:k]
        d_sub = dists[part]
        d_max = max(float(d_sub.max()), 1e-8)
        u = np.clip(d_sub / d_max, 0.0, 1.0)
        w = (1.0 - u ** 3) ** 3
        xw = xs[part]
        yw = ys[part]
        sw = w.sum()
        mx = float((w * xw).sum() / sw)
        my = float((w * yw).sum() / sw)
        sxx = float((w * (xw - mx) ** 2).sum())
        if sxx > 1e-10:
            slope = float((w * (xw - mx) * (yw - my)).sum() / sxx)
            pred[idx_g] = my + slope * (gx - mx)
        else:
            pred[idx_g] = my
    return grid, np.maximum(pred, 1e-4)


def limma_voom_fit(
    counts: np.ndarray,
    design: np.ndarray,
    coef_idx: int,
) -> dict[str, np.ndarray | float | int]:
    """Law et al. (2014) limma-voom + Smyth (2004) empirical Bayes moderated contrast."""
    counts = np.asarray(counts, dtype=np.float64)
    X = np.asarray(design, dtype=np.float64)
    n_genes, n_samples = counts.shape
    rank = int(np.linalg.matrix_rank(X))
    if rank < X.shape[1]:
        raise ValueError(f"Rank-deficient design matrix: shape={X.shape}, rank={rank}")
    df_res = n_samples - X.shape[1]
    if df_res < 2:
        raise ValueError(f"Insufficient residual df ({df_res}) for limma-voom")

    tmm = calc_tmm_norm_factors(counts)
    lib_eff = counts.sum(axis=0) * tmm
    Y = np.log2((counts + 0.5) / (lib_eff[None, :] + 1.0) * 1e6)

    # Step 1: unweighted OLS for mean-variance trend
    XtX_inv = np.linalg.inv(X.T @ X)
    B0 = Y @ (X @ XtX_inv)
    fitted_Y = B0 @ X.T
    res0 = Y - fitted_Y
    s_rsd = np.sqrt(np.maximum((res0 ** 2).sum(axis=1) / df_res, 1e-12))  # s^(1/4) = sqrt(sigma_hat) -> voom fits sqrt(sd)

    # Law et al. 2014: x-axis = average log2 count size = mean(Y) + log2(geom_mean(lib_eff + 1)) - log2(1e6)
    avg_log_count = Y.mean(axis=1) + float(np.mean(np.log2(lib_eff + 1.0))) - math.log2(1e6)
    sqrt_sd = np.sqrt(s_rsd)
    grid_x, grid_y = _lowess_1d_smooth(avg_log_count, sqrt_sd, n_grid=200, span=0.5)

    fitted_log_count = fitted_Y + np.log2(lib_eff[None, :] + 1.0) - math.log2(1e6)
    pred_sqrt_sd = np.interp(fitted_log_count, grid_x, grid_y)
    W = pred_sqrt_sd ** (-4.0)  # genes x samples precision weights

    # Step 2: vectorized per-gene weighted least squares
    # XTWX: (n_genes, p, p), XTWy: (n_genes, p)
    XTWX = np.einsum("gn,ni,nj->gij", W, X, X)
    XTWy = np.einsum("gn,ni,gn->gi", W, X, Y)
    B = np.linalg.solve(XTWX, XTWy[..., None])[..., 0]
    fitted_W = B @ X.T
    s2 = (W * (Y - fitted_W) ** 2).sum(axis=1) / df_res
    cov_B = np.linalg.inv(XTWX)
    unscaled_var = np.maximum(cov_B[:, coef_idx, coef_idx], 1e-15)

    s2_post, d0, s0 = squeeze_var(s2, df_res)
    se = np.sqrt(s2_post * unscaled_var)
    coef = B[:, coef_idx]
    t_stat = coef / se
    df_tot = float(df_res + d0)
    p_val = 2.0 * stats.t.sf(np.abs(t_stat), df=df_tot)
    q_val = bh_fdr(p_val)
    z_val = signed_z_from_t_df(t_stat, df_tot)

    return {
        "logFC": coef,
        "AveExpr": Y.mean(axis=1),
        "se": se,
        "t": t_stat,
        "z": z_val,
        "p": p_val,
        "q": q_val,
        "df_residual": int(df_res),
        "prior_df": float(d0),
        "df_total": df_tot,
        "voom_log2cpm": Y,
    }


def ebayes_moderated_ols(
    M: np.ndarray,
    design: np.ndarray,
    coef_idx: int,
    df_override: int | None = None,
) -> dict[str, np.ndarray | float | int]:
    """Vectorized empirical-Bayes moderated OLS on continuous log2 expression matrix (genes x samples)."""
    M = np.asarray(M, dtype=np.float64)
    X = np.asarray(design, dtype=np.float64)
    n_genes, n_samples = M.shape
    rank = int(np.linalg.matrix_rank(X))
    if rank < X.shape[1]:
        raise ValueError(f"Rank-deficient OLS design: shape={X.shape}, rank={rank}")
    df_ols = n_samples - X.shape[1]
    df_eff = df_ols if df_override is None else min(df_ols, int(df_override))
    if df_eff < 2:
        raise ValueError(f"Insufficient residual df ({df_eff})")

    XtX_inv = np.linalg.inv(X.T @ X)
    B = M @ (X @ XtX_inv)
    R = M - B @ X.T
    s2 = (R ** 2).sum(axis=1) / df_ols
    s2_post, d0, s0 = squeeze_var(s2, df_eff)
    se = np.sqrt(s2_post * max(float(XtX_inv[coef_idx, coef_idx]), 1e-15))
    coef = B[:, coef_idx]
    t_stat = coef / se
    df_tot = float(df_eff + d0)
    p_val = 2.0 * stats.t.sf(np.abs(t_stat), df=df_tot)
    q_val = bh_fdr(p_val)
    z_val = signed_z_from_t_df(t_stat, df_tot)
    return {
        "logFC": coef,
        "AveExpr": M.mean(axis=1),
        "se": se,
        "t": t_stat,
        "z": z_val,
        "p": p_val,
        "q": q_val,
        "df_residual": int(df_eff),
        "prior_df": float(d0),
        "df_total": df_tot,
    }


def limma_duplicate_correlation_gls(
    M: np.ndarray,
    design: np.ndarray,
    block_ids: list[str],
    coef_idx: int,
) -> dict[str, np.ndarray | float | int]:
    """Smyth et al. (2005) duplicateCorrelation block-diagonal GLS + empirical Bayes moderation.
    Properly handles mixed paired + unpaired multi-sample donors/subjects without pseudoreplication!
    """
    M = np.asarray(M, dtype=np.float64)
    X = np.asarray(design, dtype=np.float64)
    n_genes, n_samples = M.shape
    assert len(block_ids) == n_samples

    # Step 1: initial OLS residuals to estimate consensus intra-block correlation rho
    XtX_inv0 = np.linalg.inv(X.T @ X)
    B0 = M @ (X @ XtX_inv0)
    R0 = M - B0 @ X.T  # genes x samples

    blocks = defaultdict(list)
    for idx_s, b in enumerate(block_ids):
        blocks[b].append(idx_s)
    n_unique_blocks = len(blocks)

    pair_prods = []
    block_vars = []
    for idxs in blocks.values():
        if len(idxs) >= 2:
            sub = R0[:, idxs]  # genes x m
            m = len(idxs)
            # Sum of off-diagonal products per gene: ((sum_i r_i)^2 - sum_i r_i^2) / (m*(m-1))
            s_all = sub.sum(axis=1)
            s_sq = (sub ** 2).sum(axis=1)
            off = (s_all ** 2 - s_sq) / (m * (m - 1))
            pair_prods.append(off)
            block_vars.append(s_sq / m)

    if pair_prods:
        mean_cov = np.mean(np.column_stack(pair_prods), axis=1)
        tot_var = np.maximum((R0 ** 2).mean(axis=1), 1e-12)
        gene_rhos = np.clip(mean_cov / tot_var, -0.5, 0.95)
        # Consensus tanh-atanh trimmed mean (limma convention)
        atanh_r = np.arctanh(np.clip(gene_rhos, -0.9, 0.9))
        lo, hi = np.percentile(atanh_r, [15, 85])
        trimmed = atanh_r[(atanh_r >= lo) & (atanh_r <= hi)]
        rho = float(np.clip(np.tanh(np.mean(trimmed)), 0.0, 0.90))
    else:
        rho = 0.0

    # Step 2: Whiten M and X block-by-block using exact inverse square root of V_b = (1-rho)I + rho*1*1^T
    W_half = np.zeros((n_samples, n_samples), dtype=np.float64)
    for idxs in blocks.values():
        m = len(idxs)
        if m == 1 or rho <= 1e-6:
            W_half[idxs[0], idxs[0]] = 1.0
        else:
            V = (1.0 - rho) * np.eye(m) + rho * np.ones((m, m))
            eigvals, eigvecs = np.linalg.eigh(V)
            V_inv_sqrt = eigvecs @ np.diag(1.0 / np.sqrt(np.maximum(eigvals, 1e-8))) @ eigvecs.T
            W_half[np.ix_(idxs, idxs)] = V_inv_sqrt

    M_white = M @ W_half.T
    X_white = W_half @ X

    # Conservative effective degrees of freedom bounded by unique independent blocks minus fixed parameters
    df_block_bound = max(2, n_unique_blocks - X.shape[1])
    res = ebayes_moderated_ols(M_white, X_white, coef_idx=coef_idx, df_override=df_block_bound)
    res["AveExpr"] = M.mean(axis=1)
    res["consensus_rho"] = rho
    res["n_unique_blocks"] = n_unique_blocks
    return res


# ---------------------------------------------------------------------------
# 5. Global STRING Random Walk with Restart (G1 Historical Baseline)
# ---------------------------------------------------------------------------
def rwr_propagate(
    seed_weights: np.ndarray,
    left_idx: np.ndarray,
    right_idx: np.ndarray,
    edge_weights: np.ndarray,
    restart_prob: float = 0.5,
    tol: float = 1e-10,
    max_iter: int = 10000,
) -> tuple[np.ndarray, dict]:
    """Exact column-normalized symmetric adjacency RWR matching v1 G1 implementation."""
    seed = np.asarray(seed_weights, dtype=np.float64)
    if seed.max() <= 0:
        raise ValueError("All RWR seeds are zero")
    seed = seed / seed.max()
    seed = seed / seed.sum()

    n = len(seed)
    deg = np.bincount(left_idx, weights=edge_weights, minlength=n) + np.bincount(
        right_idx, weights=edge_weights, minlength=n
    )
    if np.any(deg <= 0):
        raise ValueError("Zero-degree node encountered inside induced RWR universe")

    p = seed.copy()
    for it in range(1, max_iter + 1):
        walk = np.bincount(left_idx, weights=edge_weights * p[right_idx] / deg[right_idx], minlength=n) + np.bincount(
            right_idx, weights=edge_weights * p[left_idx] / deg[left_idx], minlength=n
        )
        nxt = restart_prob * seed + (1.0 - restart_prob) * walk
        diff = float(np.abs(nxt - p).sum())
        p = nxt
        if diff < tol:
            return p, {"iterations": it, "final_l1_diff": diff, "mass_sum": float(p.sum()), "converged": True}
    raise RuntimeError(f"RWR failed to converge within {max_iter} iterations")
