#!/usr/bin/env python3
"""Render publication-quality SVG and PNG figures for EndoTarget v2 directly from
committed summary JSON/TSV artifacts, with full SHA-256 provenance receipt."""
from __future__ import annotations

import datetime as dt
import gzip
import json
import os
import pathlib
import sys

os.environ["OMP_NUM_THREADS"] = "8"
os.environ["OPENBLAS_NUM_THREADS"] = "8"
os.environ["MKL_NUM_THREADS"] = "8"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from freeze_protocol import sha256_file  # noqa: E402

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.3,
    "axes.titlesize": 10.2,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 180,
})

METHOD_COLORS = {
    "B0": "#6C7A89",  # Slate grey (naive paired DE)
    "G1": "#C0392B",  # Muted crimson (hub-biased RWR)
    "M1": "#1F77B4",  # Deep cobalt (host-decontaminated dual-contrast)
    "M2": "#16806D",  # Teal emerald (genetics + niche L-R anchored)
}


def save_both(fig: plt.Figure, out_dir: pathlib.Path, stem: str) -> list[pathlib.Path]:
    svg_p = out_dir / f"{stem}.svg"
    png_p = out_dir / f"{stem}.png"
    fig.savefig(svg_p, bbox_inches="tight")
    fig.savefig(png_p, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return [svg_p, png_p]


def render_fig1(disc_sum: dict, full_df: pd.DataFrame, cand_df: pd.DataFrame, out_dir: pathlib.Path):
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2), gridspec_kw={"width_ratios": [1.18, 1.0]})

    # Panel A: Dual-contrast scatter in GSE179640 stromal fibroblasts
    ax = axes[0]
    x = full_df["logFC_EcP_vs_EuE"].values
    y = full_df["logFC_EcP_vs_AdjP"].values
    concord = ((x > 0) & (y > 0)) | ((x < 0) & (y < 0))
    ax.scatter(x[~concord], y[~concord], s=5, c="#D5D8DC", alpha=0.45, rasterized=True, label="Host-confounded / discordant")
    ax.scatter(x[concord], y[concord], s=5, c="#AED6F1", alpha=0.45, rasterized=True, label="Dual-contrast concordant")

    b0_20 = cand_df[cand_df["in_B0_top20"].astype(str) == "True"]
    m2_20 = cand_df[cand_df["in_M2_top20"].astype(str) == "True"]
    ax.scatter(
        b0_20["logFC_EcP_vs_EuE"], b0_20["logFC_EcP_vs_AdjP"],
        s=48, facecolor=METHOD_COLORS["B0"], edgecolor="black", linewidth=0.6, zorder=4, label="B0 Top-20 (12/20 discordant!)"
    )
    ax.scatter(
        m2_20["logFC_EcP_vs_EuE"], m2_20["logFC_EcP_vs_AdjP"],
        s=54, marker="D", facecolor=METHOD_COLORS["M2"], edgecolor="black", linewidth=0.6, zorder=5, label="M2 Top-20 (20/20 concordant)"
    )

    annot_genes = {"LCN2", "GPC3", "S100A10", "SFRP2", "PITX1", "HMGA2", "CXCL8", "KCNK3", "MMP1", "SPHK1", "CDKN2B", "CXCR4"}
    for _, row in pd.concat([b0_20, m2_20]).drop_duplicates(subset=["ensembl_gene_id"]).iterrows():
        sym = row["symbol"]
        if sym in annot_genes:
            dx = 0.12 if row["logFC_EcP_vs_EuE"] >= 0 else -0.15
            dy = 0.13 if row["logFC_EcP_vs_AdjP"] >= 0 else -0.24
            ax.text(row["logFC_EcP_vs_EuE"] + dx, row["logFC_EcP_vs_AdjP"] + dy, sym, fontsize=7.6, fontweight="bold")

    ax.axhline(0, color="#444444", lw=0.9, ls="--")
    ax.axvline(0, color="#444444", lw=0.9, ls="--")
    ax.set_xlabel("Discovery Paired log2FC: Ectopic Peritoneal Lesion vs Eutopic (n=8 pairs)")
    ax.set_ylabel("Discovery Paired log2FC: Ectopic Lesion vs Adjacent Peritoneum (n=6 triplets)")
    ax.set_title("A. Dual-Contrast Geometry Separates Lesion Remodeling from Peritoneal Background")
    ax.legend(loc="lower right", frameon=True, fontsize=7.8)

    # Panel B: Top-20 shortlist diagnostic profile across B0, G1, M1, M2
    ax2 = axes[1]
    methods = ["B0", "G1", "M1", "M2"]
    diag = disc_sum["method_shortlist_diagnostics"]
    categories = [
        ("Dual-contrast\nconcordant (%)", [diag[m]["top20"]["dual_concordant_count"] * 5 for m in methods]),
        ("Mesothelial\nambient >4x (%)", [diag[m]["top20"]["ambient_meso_enriched_gt2_count"] * 5 for m in methods]),
        ("Human genetics\nsupported (%)", [diag[m]["top20"]["human_genetics_supported_count"] * 5 for m in methods]),
        ("Active niche\nL-R axis (%)", [diag[m]["top20"]["active_niche_lr_count"] * 5 for m in methods]),
        ("SM / AB\ntractable (%)", [diag[m]["top20"]["tractable_sm_or_ab_count"] * 5 for m in methods]),
    ]
    x_pos = np.arange(len(categories))
    width = 0.19
    for idx_m, m in enumerate(methods):
        vals = [c[1][idx_m] for c in categories]
        bars = ax2.bar(
            x_pos + (idx_m - 1.5) * width, vals, width=width,
            color=METHOD_COLORS[m], edgecolor="black", linewidth=0.5,
            label=f"{m} (med STRING deg={diag[m]['top20']['median_string_degree']:.1f})",
        )
        for b in bars:
            h = b.get_height()
            if h > 0:
                ax2.text(b.get_x() + b.get_width() / 2.0, h + 1.5, f"{int(round(h))}", ha="center", va="bottom", fontsize=7.1)

    ax2.set_xticks(x_pos)
    ax2.set_xticklabels([c[0] for c in categories], fontsize=8.1)
    ax2.set_ylim(0, 114)
    ax2.set_ylabel("Percentage of Frozen Top-20 Candidates (%)")
    ax2.set_title("B. Pre-Seal Discovery Audit of Frozen Top-20 Candidate Shortlists")
    ax2.legend(loc="upper right", frameon=True, fontsize=7.8)

    fig.tight_layout()
    return save_both(fig, out_dir, "fig1_discovery_dual_contrast_and_shortlist_audit")


def render_fig2(eval_sum: dict, out_dir: pathlib.Path):
    sc_res = eval_sum["GSE213216_external_scrna"]["primary_stromal_EcP_vs_EuE"]
    sens = eval_sum["GSE213216_external_scrna"]["sensitivities"]
    sm = sc_res["shortlist_metrics"]
    boot = sc_res["donor_cluster_bootstrap"]["comparisons"]

    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.7))
    methods = ["B0", "G1", "M1", "M2"]

    # Panel A: GSE213216 EcP vs EuE (uncontrolled peritoneal-vs-uterine scRNA) vs Endometrioma EcO vs EuE
    ax = axes[0]
    x = np.arange(len(methods))
    w = 0.25
    rep_ecp = [sm[m]["top20"]["replicated_q05_rate_on_k"] * 100 for m in methods]
    dir_ecp = [sm[m]["top20"]["same_direction_rate_evaluable"] * 100 for m in methods]
    dir_eco = [sens["Endometrioma_EcO_vs_EuE"]["shortlists"][m]["top20"]["same_direction_rate_evaluable"] * 100 for m in methods]

    ax.bar(x - w, rep_ecp, width=w, color=[METHOD_COLORS[m] for m in methods], edgecolor="black", lw=0.6, label="EcP vs EuE q<0.05 (%)")
    ax.bar(x, dir_ecp, width=w, color="#AED6F1", edgecolor="black", lw=0.6, label="EcP vs EuE Same Sign (%)")
    ax.bar(x + w, dir_eco, width=w, color="#D7BDE2", edgecolor="black", lw=0.6, label="EcO vs EuE Same Sign (%)")
    ax.axhline(50, color="#666666", ls=":", lw=0.9)

    for i in range(len(methods)):
        ax.text(x[i] - w, rep_ecp[i] + 1.5, f"{rep_ecp[i]:.0f}%", ha="center", fontsize=7.3, fontweight="bold")
        ax.text(x[i], dir_ecp[i] + 1.5, f"{dir_ecp[i]:.0f}%", ha="center", fontsize=7.1)
        ax.text(x[i] + w, dir_eco[i] + 1.5, f"{dir_eco[i]:.0f}%", ha="center", fontsize=7.1)

    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontweight="bold")
    ax.set_ylim(0, 112)
    ax.set_ylabel("External scRNA Rate (%)")
    ax.set_title("A. GSE213216 scRNA: Peritoneal & Ovarian Lesions vs EuE")
    ax.legend(loc="upper right", fontsize=7.3, frameon=True)

    # Panel B: Why B0 'replicates' in EcP vs EuE yet collapses in Host Peritoneum (Top-20 B0 gene trajectories)
    ax2 = axes[1]
    top_audit = pd.read_csv(ROOT / "reports/top20_external_audit_v2.tsv", sep="\t")
    b0_sub = top_audit[(top_audit["method"] == "B0") & (top_audit["GSE141549_PeL_vs_PP_HostCtrl_logFC"].notna())].copy()
    b0_sub["GSE141549_PeL_vs_PE_logFC"] = b0_sub["GSE141549_PeL_vs_PE_logFC"].astype(float)
    b0_sub["GSE141549_PeL_vs_PP_HostCtrl_logFC"] = b0_sub["GSE141549_PeL_vs_PP_HostCtrl_logFC"].astype(float)
    sc_x = b0_sub["GSE213216_logFC"].values * b0_sub["discovery_sign"].values
    pp_y = b0_sub["GSE141549_PeL_vs_PP_HostCtrl_logFC"].values * b0_sub["discovery_sign"].values
    ax2.scatter(sc_x, pp_y, s=52, c=METHOD_COLORS["B0"], edgecolor="black", lw=0.7)
    for _, r in b0_sub.iterrows():
        sx = float(r["GSE213216_logFC"]) * int(r["discovery_sign"])
        py = float(r["GSE141549_PeL_vs_PP_HostCtrl_logFC"]) * int(r["discovery_sign"])
        if abs(sx) > 2.8 or abs(py) > 0.65:
            dx_lbl = -1.55 if r["symbol"] == "S100A10" else 0.16
            dy_lbl = 0.05 if r["symbol"] == "S100A10" else (0.04 if py >= 0 else -0.05)
            ax2.text(sx + dx_lbl, py + dy_lbl, r["symbol"], fontsize=7.4, fontweight="bold")
    ax2.axhline(0, color="#C0392B", ls="--", lw=1.1)
    ax2.axvline(0, color="#555555", ls="--", lw=0.9)
    ax2.set_xlabel("GSE213216 scRNA Directional log2FC (EcP vs EuE)")
    ax2.set_ylabel("GSE141549 Host-Controlled Directional log2FC (PeL vs PP)")
    ax2.set_title("B. B0's scRNA Hits Invert Against Host Peritoneum (PP)")

    # Panel C: Donor-Cluster Bootstrap Differences in GSE213216 EcP vs EuE (2,000 draws)
    ax3 = axes[2]
    comps = ["M1_minus_B0", "M1_minus_G1", "M2_minus_B0", "M2_minus_G1"]
    labels_c = ["M1 − B0 (Top-20)", "M1 − G1 (Top-20)", "M2 − B0 (Top-20)", "M2 − G1 (Top-20)",
                "M1 − B0 (Top-50)", "M1 − G1 (Top-50)", "M2 − B0 (Top-50)", "M2 − G1 (Top-50)"]
    y_pos = np.arange(len(labels_c))[::-1]
    means, lows, highs = [], [], []
    for k_tag in ("top20", "top50"):
        for c in comps:
            d = boot[c][k_tag]["rep_rate_on_k"]
            means.append(d["boot_mean"] * 100)
            lows.append(d["ci_95"][0] * 100)
            highs.append(d["ci_95"][1] * 100)

    ax3.axvline(0, color="#555555", ls="--", lw=1.0)
    for idx_y, (yp, mn, lo, hi) in enumerate(zip(y_pos, means, lows, highs)):
        col = METHOD_COLORS["M1"] if "M1" in labels_c[idx_y] else METHOD_COLORS["M2"]
        ax3.errorbar(
            mn, yp, xerr=[[mn - lo], [hi - mn]],
            fmt="o", color=col, ecolor=col, elinewidth=2.0, capsize=4, markersize=5.5
        )
        ax3.text(hi + 1.5, yp, f"{mn:+.1f} [{lo:+.1f}, {hi:+.1f}]", va="center", fontsize=7.3)

    ax3.set_yticks(y_pos)
    ax3.set_yticklabels(labels_c, fontsize=8.1)
    ax3.set_xlim(min(lows) - 8, max(highs) + 46)
    ax3.set_xlabel("Δ Replicated Rate in EcP vs EuE (pp, 95% CI)")
    ax3.set_title("C. GSE213216 Donor-Cluster Bootstrap (n=17 Patients)")

    fig.tight_layout()
    return save_both(fig, out_dir, "fig2_external_scrna_replication_GSE213216")


def render_fig3(eval_sum: dict, out_dir: pathlib.Path):
    bulk = eval_sum["GSE141549_external_bulk"]
    j20 = {m: bulk["Joint_Dual_Contrast_and_Host_False_Positive_Audit"][m]["top20"] for m in ("B0", "G1", "M1", "M2")}
    p_pp = bulk["Paired_Peritoneal_Lesion_vs_PP_HostControlled_n20"]["shortlist_metrics"]
    prog_host = bulk["Program_Transfer_All_Peritoneal_vs_Host"]
    sub = bulk["Multi_Subtype_and_Control_Contrasts"]
    methods = ["B0", "G1", "M1", "M2"]

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.9))

    # Panel A: Paired Origin Contrast (PeL vs PE, n=41) Collapse When Tested Against Host Peritoneum (PeL vs PP, n=20)
    ax = axes[0]
    x = np.arange(len(methods))
    w = 0.25
    dir_pe = [bulk["Paired_Peritoneal_Lesion_vs_PE_n41"]["shortlist_metrics"][m]["top20"]["same_direction_rate_evaluable"] * 100 for m in methods]
    dir_pp = [p_pp[m]["top20"]["same_direction_rate_evaluable"] * 100 for m in methods]
    dir_die_pp = [sub["Paired_DIE_vs_PP_n20"]["shortlists"][m]["top20"]["same_direction_rate_evaluable"] * 100 for m in methods]

    ax.bar(x - w, dir_pe, width=w, color="#BDC3C7", edgecolor="black", lw=0.6, label="Paired PeL vs PE (n=41)")
    ax.bar(x, dir_pp, width=w, color=[METHOD_COLORS[m] for m in methods], edgecolor="black", lw=0.7, label="Paired PeL vs PP Host-Ctrl (n=20)")
    ax.bar(x + w, dir_die_pp, width=w, color=[METHOD_COLORS[m] for m in methods], alpha=0.65, hatch="//", edgecolor="black", lw=0.6, label="Paired DIE vs PP Host-Ctrl (n=20)")
    ax.axhline(50, color="#555555", ls=":", lw=1.0)

    for i in range(len(methods)):
        ax.text(x[i] - w, dir_pe[i] + 1.5, f"{dir_pe[i]:.0f}%", ha="center", fontsize=7.0)
        ax.text(x[i], dir_pp[i] + 1.5, f"{dir_pp[i]:.0f}%", ha="center", fontsize=7.4, fontweight="bold")
        ax.text(x[i] + w, dir_die_pp[i] + 1.5, f"{dir_die_pp[i]:.0f}%", ha="center", fontsize=7.1)

    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontweight="bold")
    ax.set_ylim(0, 114)
    ax.set_ylabel("Direction Concordance (% Evaluable, Top-20)")
    ax.set_title("A. Host-Peritoneum Control Reverses Ranking Hierarchy")
    ax.legend(loc="upper right", fontsize=7.3, frameon=True)

    # Panel B: Host-Controlled Program Transfer (All Peritoneal Lesion n=54 vs Host Peritoneum PP+CP n=61)
    ax2 = axes[1]
    thetas = [prog_host[m]["top20"]["theta"] for m in methods]
    lo = [prog_host[m]["top20"]["ci_95"][0] for m in methods]
    hi = [prog_host[m]["top20"]["ci_95"][1] for m in methods]
    aucs = [prog_host[m]["top20"]["auc"] for m in methods]
    pvals = [prog_host[m]["top20"]["p_value"] for m in methods]

    ax2.axhline(0, color="#555555", ls="--", lw=1.0)
    for i, m in enumerate(methods):
        ax2.errorbar(
            i, thetas[i], yerr=[[thetas[i] - lo[i]], [hi[i] - thetas[i]]],
            fmt="D", color=METHOD_COLORS[m], ecolor=METHOD_COLORS[m], elinewidth=2.2, capsize=5, markersize=7.5
        )
        y_txt = hi[i] + 0.014 if thetas[i] >= 0 else lo[i] - 0.052
        ax2.text(i, y_txt, f"AUC={aucs[i]:.2f}\np={pvals[i]:.3g}", ha="center", fontsize=7.6, fontweight="bold")

    ax2.set_xticks(x)
    ax2.set_xticklabels(methods, fontweight="bold")
    ax2.set_ylim(-0.23, 0.25)
    ax2.set_ylabel("Program Shift θ vs Host Peritoneum (PP+CP)")
    ax2.set_title("B. Host-Controlled Program Transfer (n=54 vs 61 Subjects)")

    # Panel C: Multi-Subtype Mean Directional Z Across Host-Controlled Bulk Contrasts (Top-20)
    ax3 = axes[2]
    z_pel_pp = [p_pp[m]["top20"]["mean_directional_z"] for m in methods]
    z_die_all = [sub["All_DIE_vs_PP_CP"]["shortlists"][m]["top20"]["mean_directional_z"] for m in methods]
    z_oma_all = [sub["All_OMA_vs_PP_CP"]["shortlists"][m]["top20"]["mean_directional_z"] for m in methods]

    ax3.axhline(0, color="#555555", ls="--", lw=0.9)
    ax3.bar(x - w, z_pel_pp, width=w, color="#2E86C1", edgecolor="black", lw=0.6, label="Paired PeL vs PP (n=20)")
    ax3.bar(x, z_die_all, width=w, color="#17A589", edgecolor="black", lw=0.6, label="DIE vs PP+CP (n=57/61)")
    ax3.bar(x + w, z_oma_all, width=w, color="#8E44AD", edgecolor="black", lw=0.6, label="OMA vs PP+CP (n=25/61)")

    for i in range(len(methods)):
        for j_bar, (offset, val) in enumerate(((-w, z_pel_pp[i]), (0.0, z_die_all[i]), (w, z_oma_all[i]))):
            extra = 0.14 if (i == 1 and j_bar == 2) else 0.0
            yt = val + 0.07 + extra if val >= 0 else val - 0.23
            ax3.text(x[i] + offset, yt, f"{val:+.2f}", ha="center", fontsize=6.7, fontweight="bold")

    ax3.set_xticks(x)
    ax3.set_xticklabels(methods, fontweight="bold")
    ax3.set_ylim(-2.7, 2.9)
    ax3.set_ylabel("Mean Directional External Z-Statistic (Top-20)")
    ax3.set_title("C. Multi-Subtype Host-Controlled Directional Z (GSE141549)")
    ax3.legend(loc="lower right", fontsize=7.3, frameon=True)

    fig.tight_layout()
    return save_both(fig, out_dir, "fig3_bulk_double_difference_and_subtypes_GSE141549")


def render_fig4(eval_sum: dict, panel_path: pathlib.Path, out_dir: pathlib.Path):
    gen_enr = eval_sum["human_genetics_enrichment"]
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 5.6), gridspec_kw={"width_ratios": [1.02, 1.35]})

    # Panel A: Compartment-resolved genetic enrichment (Mean Rank Percentile & Mann-Whitney p)
    ax = axes[0]
    comp_keys = [
        ("GSE141549_Paired_PeL_vs_PP_HostControlled_absZ", "GSE141549 Paired PeL vs PP"),
        ("GSE213216_stromal_EcP_vs_EuE_absZ", "GSE213216 scRNA Stroma"),
        ("GSE179640_perivascular_mural_absZ", "GSE179640 Perivascular Mural"),
        ("GSE179640_endothelial_absZ", "GSE179640 Endothelial"),
        ("GSE179640_epithelial_absZ", "GSE179640 Epithelial"),
        ("GSE179640_stromal_fibroblast_absZ", "GSE179640 Stroma (Naive)"),
        ("GSE179640_myeloid_macrophage_absZ", "GSE179640 Myeloid / Macrophage"),
    ]
    present_ck = [(k, lbl) for k, lbl in comp_keys if k in gen_enr]
    y_pos = np.arange(len(present_ck))[::-1]
    pct_ot = [gen_enr[k]["OpenTargets_Genetic_ge0.5"]["mean_rank_percentile"] for k, _ in present_ck]
    p_ot = [gen_enr[k]["OpenTargets_Genetic_ge0.5"]["mannwhitney_p_one_sided"] for k, _ in present_ck]
    pct_rah = [gen_enr[k]["Rahmioglu2023_CredibleNearest"]["mean_rank_percentile"] for k, _ in present_ck]
    p_rah = [gen_enr[k]["Rahmioglu2023_CredibleNearest"]["mannwhitney_p_one_sided"] for k, _ in present_ck]

    ax.axvline(50.0, color="#666666", ls="--", lw=1.0, label="Null expectation (50th pct)")
    ax.barh(y_pos + 0.17, pct_ot, height=0.32, color="#16806D", edgecolor="black", lw=0.5, label="Open Targets Genetic >=0.5")
    ax.barh(y_pos - 0.17, pct_rah, height=0.32, color="#2E86C1", edgecolor="black", lw=0.5, label="Rahmioglu 2023 Credible Set")
    for yp, val, pv in zip(y_pos, pct_ot, p_ot):
        sig = "**" if pv < 0.01 else ("*" if pv < 0.05 else "")
        ax.text(val + 0.5, yp + 0.17, f"{val:.1f}% (p={pv:.2g}){sig}", va="center", fontsize=7.2)

    ax.set_yticks(y_pos)
    ax.set_yticklabels([lbl for _, lbl in present_ck], fontsize=8.2)
    ax.set_xlim(38, max(max(pct_ot), max(pct_rah)) + 16)
    ax.set_xlabel("Mean Lesion Remodeling Rank Percentile (%)")
    ax.set_title("A. Compartment & Cohort Human Genetic Enrichment")
    ax.legend(loc="lower right", fontsize=7.5, frameon=True)

    # Panel B: Multi-Evidence Dossier & Candidate Evidence Matrix (12 targets)
    ax2 = axes[1]
    ddf = pd.read_csv(panel_path, sep="\t")
    dsign = ddf["discovery_sign"].values
    genes_lbl = [
        f"{'[Tier 1] ' if str(r['evidence_category']).startswith('Tier_1') else ''}{r['symbol']} (M2 #{r['rank_M2']})"
        for _, r in ddf.iterrows()
    ]
    # Signed directional evidence scaled to [-1, +1] for expression contrasts, [0, 1] for priors/tractability
    mat = np.column_stack([
        np.clip((ddf["disc_logFC_EcP_vs_EuE"].values * dsign) / 3.0, -1, 1),
        np.clip((ddf["disc_logFC_EcP_vs_AdjP"].values * dsign) / 2.0, -1, 1),
        np.clip((ddf["ext_sc_GSE213216_logFC"].values * dsign) / 2.5, -1, 1),
        np.clip((ddf["ext_bulk_PeL_vs_PP_HostCtrl_logFC"].values * dsign) / 1.0, -1, 1),
        np.clip(ddf["genetic_prior_score"].values, 0, 1),
        np.clip(ddf["niche_lr_score"].values, 0, 1),
        ((ddf["tractable_sm"].astype(str) == "True") | (ddf["tractable_ab"].astype(str) == "True")).astype(float),
    ])
    im = ax2.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-1.0, vmax=1.0)
    col_names = [
        "Disc\nEcP-EuE", "Disc\nEcP-AdjP", "Ext scRNA\nGSE213216",
        "Host Bulk\nPeL-PP", "Human\nGenetics", "Niche\nL-R", "SM/AB\nTract"
    ]
    ax2.set_xticks(range(len(col_names)))
    ax2.set_xticklabels(col_names, fontsize=7.9)
    ax2.set_yticks(range(len(genes_lbl)))
    ax2.set_yticklabels(genes_lbl, fontsize=8.2, fontweight="bold")
    for r_i in range(mat.shape[0]):
        for c_i in range(mat.shape[1]):
            val = mat[r_i, c_i]
            ax2.text(
                c_i, r_i, f"{val:+.2f}" if c_i < 4 else f"{val:.2f}",
                ha="center", va="center", fontsize=6.9,
                color="white" if abs(val) > 0.58 else "#1B2631",
            )
    plt.colorbar(im, ax=ax2, fraction=0.035, pad=0.02, label="Directional Concordance [-1, +1] / Prior [0, 1]")
    ax2.set_title("B. Multi-Evidence Target Dossier & Candidate Panel")

    fig.tight_layout()
    return save_both(fig, out_dir, "fig4_genetics_enrichment_and_target_dossiers")


def main():
    disc_sum_p = ROOT / "reports/discovery/discovery_summary_v2.json"
    eval_sum_p = ROOT / "reports/locked_evaluation_summary_v2.json"
    full_gz_p = ROOT / "reports/discovery/GSE179640_stromal_rankings_full.tsv.gz"
    cand_p = ROOT / "manifests/frozen_candidates_v2.tsv"
    panel_p = ROOT / "reports/target_dossiers_extended_panel.tsv"

    disc_sum = json.loads(disc_sum_p.read_text(encoding="utf-8"))
    eval_sum = json.loads(eval_sum_p.read_text(encoding="utf-8"))
    with gzip.open(full_gz_p, "rt", encoding="utf-8") as f:
        full_df = pd.read_csv(f, sep="\t")
    cand_df = pd.read_csv(cand_p, sep="\t")

    out_dir = ROOT / "reports/figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    rendered = []
    rendered.extend(render_fig1(disc_sum, full_df, cand_df, out_dir))
    rendered.extend(render_fig2(eval_sum, out_dir))
    rendered.extend(render_fig3(eval_sum, out_dir))
    rendered.extend(render_fig4(eval_sum, panel_p, out_dir))

    receipt = {
        "schema_version": 2,
        "phase": "STAGE3_FIGURES_RENDER",
        "rendered_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_artifacts": {
            str(p.relative_to(ROOT)): sha256_file(p)
            for p in (disc_sum_p, eval_sum_p, full_gz_p, cand_p, panel_p) if p.exists()
        },
        "script_sha256": sha256_file(pathlib.Path(__file__)),
        "figures": {
            str(p.relative_to(ROOT)): {"sha256": sha256_file(p), "bytes": p.stat().st_size}
            for p in sorted(rendered)
        },
    }
    rec_p = ROOT / "receipts/figures_v2.json"
    rec_p.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "FIGURES_RENDERED", "count": len(rendered), "receipt_sha256": sha256_file(rec_p)}, indent=2))


if __name__ == "__main__":
    main()
