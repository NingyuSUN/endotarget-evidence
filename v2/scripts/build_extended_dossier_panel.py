#!/usr/bin/env python3
"""Build descriptive 12-gene extended multi-evidence target dossier context table
(Tier 1 RUNX1 + top host-controlled replicated M2 niche targets + GWAS/coding loci replicating in external cohorts)
with GTEx v10 Uterus/Fibroblast/Ovary expression and eQTL annotations via gtex_cli.py.
Does NOT modify any locked evaluation table or receipt."""
from __future__ import annotations

import csv
import gzip
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from evaluate_locked import query_gtex_for_dossier_genes  # noqa: E402
from m0_audit_and_priors import write_tsv  # noqa: E402


def main():
    with gzip.open(ROOT / "reports/discovery/GSE179640_stromal_rankings_full.tsv.gz", "rt", encoding="utf-8") as f:
        disc = {r["symbol"]: r for r in csv.DictReader(f, delimiter="\t") if r["symbol"]}
    with (ROOT / "reports/external_scrna/GSE213216_stromal_DE_EcP_vs_EuE.tsv").open("r", encoding="utf-8") as f:
        sc213 = {r["ensembl_gene_id"]: r for r in csv.DictReader(f, delimiter="\t")}
    with (ROOT / "reports/external_bulk/GSE141549_Paired_Peritoneal_Lesion_vs_PP_HostControlled.tsv").open("r", encoding="utf-8") as f:
        b_pp = {r["ensembl_gene_id"]: r for r in csv.DictReader(f, delimiter="\t")}
    with (ROOT / "reports/external_bulk/GSE141549_Paired_Peritoneal_Lesion_vs_PE.tsv").open("r", encoding="utf-8") as f:
        b_pe = {r["ensembl_gene_id"]: r for r in csv.DictReader(f, delimiter="\t")}
    with (ROOT / "manifests/human_genetics_prior_v2.tsv").open("r", encoding="utf-8") as f:
        gen = {r["ensembl_gene_id"]: r for r in csv.DictReader(f, delimiter="\t")}

    panel_spec = [
        ("RUNX1", "Tier_1_Locked_MultiModal"),
        ("CXCR4", "M2_Top20_HostBulk_q05_NicheLR"),
        ("IL1B", "M2_Top20_HostBulk_Borderline_Genetics_NicheLR"),
        ("CXCL8", "M2_Top20_HostBulk_Borderline_NicheLR"),
        ("VMP1", "M2_Top50_HostBulk_q05"),
        ("GNG4", "M2_Top50_HostBulk_q05"),
        ("TPM4", "M2_Top50_HostBulk_q05"),
        ("SLC4A7", "M2_Top50_HostBulk_q05"),
        ("PDLIM5", "Exploratory_GWAS_CredSet_scRNA_and_HostBulk_q05"),
        ("F5", "Exploratory_CodingCredSet_rs6025_HostBulk_q05"),
        ("SMAD3", "Exploratory_GWAS_15q22_M1_and_HostBulk_q05"),
        ("ESR1", "M2_Top50_Canonical_GWAS_6q25_HormoneReceptor"),
    ]

    gtex_map = query_gtex_for_dossier_genes([s for s, _ in panel_spec])
    rows = []
    for sym, category in panel_spec:
        d = disc[sym]
        gid = d["ensembl_gene_id"]
        gr = gen.get(gid, {})
        s_rec = sc213.get(gid, {})
        pp_rec = b_pp.get(gid, {})
        pe_rec = b_pe.get(gid, {})
        gt = gtex_map.get(sym, {})
        rows.append({
            "evidence_category": category,
            "symbol": sym,
            "ensembl_gene_id": gid,
            "discovery_sign": int(d["discovery_sign"]),
            "rank_M2": int(d["rank_M2"]),
            "rank_M1": int(d["rank_M1"]),
            "rank_B0": int(d["rank_B0"]),
            "disc_logFC_EcP_vs_EuE": float(d["logFC_EcP_vs_EuE"]),
            "disc_logFC_EcP_vs_AdjP": float(d["logFC_EcP_vs_AdjP"]),
            "disc_m1_core": float(d["m1_core"]),
            "ext_sc_GSE213216_logFC": float(s_rec.get("logFC_EcP_vs_EuE", 0.0)),
            "ext_sc_GSE213216_q": float(s_rec.get("q_bh", 1.0)),
            "ext_bulk_PeL_vs_PP_HostCtrl_logFC": float(pp_rec.get("logFC", 0.0)),
            "ext_bulk_PeL_vs_PP_HostCtrl_q": float(pp_rec.get("q_bh", 1.0)),
            "ext_bulk_PeL_vs_PE_logFC": float(pe_rec.get("logFC", 0.0)),
            "ext_bulk_PeL_vs_PE_q": float(pe_rec.get("q_bh", 1.0)),
            "genetic_prior_score": float(d["genetic_prior_score"]),
            "ot_genetic_association": float(gr.get("ot_genetic_association", 0.0)),
            "niche_lr_score": float(d["niche_lr_score"]),
            "best_lr_pair": d["best_lr_pair"] or "-",
            "tractable_sm": gr.get("tractable_sm_any", "False"),
            "tractable_ab": gr.get("tractable_ab_any", "False"),
            "gtex_tpm_uterus": gt.get("tpm_Uterus", "NA"),
            "gtex_tpm_fibroblasts": gt.get("tpm_Fibroblasts", "NA"),
            "gtex_tpm_ovary": gt.get("tpm_Ovary", "NA"),
            "gtex_eqtls_uterus_fibro_ovary": (
                f"{gt.get('eqtls_Uterus', 0)}/{gt.get('eqtls_Fibroblasts', 0)}/{gt.get('eqtls_Ovary', 0)}"
                if gt else "NA"
            ),
        })

    out_p = ROOT / "reports/target_dossiers_extended_panel.tsv"
    write_tsv(out_p, rows, list(rows[0].keys()))
    print(f"Wrote {len(rows)} rows to {out_p}")


if __name__ == "__main__":
    main()
