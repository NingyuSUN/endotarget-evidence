#!/usr/bin/env python3
"""Stage 1 (M0): Metadata-only cohort audit, human genetics & tractability prior table,
de novo canonical lineage markers, curated niche L-R pairs, and STRING v12.5 direct gene network.

STRICT ZERO-LEAKAGE GUARANTEE:
- Reads ONLY GEO SOFT metadata headers (explicitly skipping any !sample_table_begin..end blocks),
  GEO supplementary file listings, GSE141549_Sample_link.xlsx (sample-to-subject link table),
  Open Targets MONDO_0005133 (EFO:0001065) target/drug/GWAS JSONs, Rahmioglu et al. 2023 Nat Genet
  supplementary tables (41588_2023_1323_MOESM4_ESM.xlsx), HGNC complete gene annotation, and
  STRING v12.5 human direct experimental/database interactions.
- Zero cohort gene-expression values are downloaded or read in Stage 1.
"""
from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import os
import pathlib
import re
from collections import Counter, defaultdict

import openpyxl

os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")

ROOT = pathlib.Path(__file__).resolve().parents[1]
ENSG_RE = re.compile(r"^ENSG[0-9]{11}$")


def sha256_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: pathlib.Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def load_hgnc(path: pathlib.Path):
    approved = {}
    prev_map = {}
    alias_map = {}
    prev_amb = set()
    alias_amb = set()
    entrez_to_ensg = {}
    ensg_to_symbol = {}
    approved_multi = defaultdict(list)

    with path.open("r", encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("status") != "Approved":
                continue
            sym = r["symbol"].strip()
            eg = r["ensembl_gene_id"].strip()
            approved_multi[sym].append(r)
            if not ENSG_RE.fullmatch(eg):
                continue
            approved[sym] = eg
            ensg_to_symbol[eg] = sym
            entrez = r.get("entrez_id", "").strip()
            if entrez:
                entrez_to_ensg[entrez] = eg
            for col, store, amb in (("prev_symbol", prev_map, prev_amb), ("alias_symbol", alias_map, alias_amb)):
                raw = r.get(col, "")
                for tok in filter(None, (x.strip().strip('"') for x in raw.split("|"))):
                    if tok in store and store[tok] != eg:
                        amb.add(tok)
                    store[tok] = eg

    prev_clean = {k: v for k, v in prev_map.items() if k not in prev_amb and k not in approved}
    alias_clean = {k: v for k, v in alias_map.items() if k not in alias_amb and k not in approved and k not in prev_clean}
    return {
        "approved": approved,
        "prev": prev_clean,
        "alias": alias_clean,
        "entrez_to_ensg": entrez_to_ensg,
        "ensg_to_symbol": ensg_to_symbol,
        "approved_multi": approved_multi,
    }


def resolve_symbol(sym: str, hgnc: dict) -> tuple[str | None, str]:
    s = sym.strip()
    if s in hgnc["approved"]:
        return hgnc["approved"][s], "approved"
    if s in hgnc["prev"]:
        return hgnc["prev"][s], "previous"
    if s in hgnc["alias"]:
        return hgnc["alias"][s], "alias"
    return None, "unmapped"


def parse_soft_metadata_only(path: pathlib.Path):
    """Parse GEO family SOFT file headers while strictly skipping any data tables."""
    series = defaultdict(list)
    samples = []
    cur = None
    kind = None
    in_data_table = False
    skipped_table_lines = 0

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for raw_line in f:
            line = raw_line.rstrip("\r\n")
            if line.startswith(("!sample_table_begin", "!platform_table_begin")):
                in_data_table = True
                continue
            if line.startswith(("!sample_table_end", "!platform_table_end")):
                in_data_table = False
                continue
            if in_data_table:
                skipped_table_lines += 1
                continue
            if line.startswith("^SERIES"):
                kind = "series"
                cur = None
                continue
            if line.startswith("^PLATFORM"):
                kind = "platform"
                cur = None
                continue
            if line.startswith("^SAMPLE"):
                kind = "sample"
                gsm = line.split("=", 1)[1].strip()
                cur = {"gsm": gsm, "fields": defaultdict(list), "chars": {}}
                samples.append(cur)
                continue
            if line.startswith("!") and "=" in line:
                k, v = line[1:].split("=", 1)
                k, v = k.strip(), v.strip()
                if kind == "series":
                    series[k].append(v)
                elif kind == "sample" and cur is not None:
                    cur["fields"][k].append(v)
                    if k == "Sample_characteristics_ch1" and ":" in v:
                        ck, cv = v.split(":", 1)
                        cur["chars"][ck.strip().lower()] = cv.strip()
    return series, samples, skipped_table_lines


def parse_filelist(path: pathlib.Path) -> dict[str, int]:
    sizes = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) >= 5 and parts[0] == "File":
            sizes[parts[1].strip()] = int(parts[3].strip())
    return sizes


def audit_gse179640(meta_dir: pathlib.Path):
    _, samples, skipped = parse_soft_metadata_only(meta_dir / "GSE179640_family.soft.gz")
    sizes = parse_filelist(meta_dir / "filelist.txt")

    sample_rows = []
    by_donor = defaultdict(lambda: defaultdict(list))

    for s in sorted(samples, key=lambda x: x["gsm"]):
        f = s["fields"]
        c = s["chars"]
        supp_urls = [u for k, vals in f.items() if k.startswith("Sample_supplementary_file") for u in vals if u != "NONE"]
        supp_names = [u.split("/")[-1] for u in supp_urls]
        method = c.get("method", "")
        lib_type = c.get("library type", "")
        tissue = c.get("tissue", "")
        loc = c.get("sample location", "")
        donor = f.get("Sample_source_name_ch1", [""])[0]
        cond = c.get("condition", "")

        is_tissue_scrna = (
            method == "scRNA-seq"
            and lib_type == "Gene Expression"
            and tissue != "Patient-Derived Organoid"
        )
        site_code = "EXCLUDED"
        if is_tissue_scrna:
            if loc == "Eutopic" and cond == "Endometriosis":
                site_code = "EuE"
            elif loc == "Eutopic" and cond == "Control":
                site_code = "Ctrl_EuE"
            elif loc == "Ectopic" and tissue == "Peritoneum":
                site_code = "EcP"
            elif loc == "Ectopic Adjacent" and tissue == "Peritoneum":
                site_code = "AdjP"
            elif loc == "Ectopic Ovary":
                site_code = "EcO"
            by_donor[donor][site_code].append(s["gsm"])

        bytes_total = sum(sizes.get(fn, 0) for fn in supp_names)
        sample_rows.append({
            "gsm": s["gsm"],
            "donor_id": donor,
            "condition": cond,
            "method": method,
            "library_type": lib_type,
            "sample_location": loc,
            "tissue": tissue,
            "site_code": site_code,
            "included_discovery_scrna": str(is_tissue_scrna),
            "supplementary_file": ";".join(supp_names),
            "supplementary_url": ";".join(supp_urls),
            "bytes": bytes_total,
        })

    donor_rows = []
    for donor in sorted(by_donor.keys()):
        sites = by_donor[donor]
        has_eue = len(sites["EuE"]) > 0
        has_ecp = len(sites["EcP"]) > 0
        has_adjp = len(sites["AdjP"]) > 0
        has_eco = len(sites["EcO"]) > 0
        donor_rows.append({
            "donor_id": donor,
            "group": "Control" if donor.startswith("C") else "Endometriosis",
            "Ctrl_EuE_gsms": ";".join(sites["Ctrl_EuE"]),
            "EuE_gsms": ";".join(sites["EuE"]),
            "EcP_gsms": ";".join(sites["EcP"]),
            "AdjP_gsms": ";".join(sites["AdjP"]),
            "EcO_gsms": ";".join(sites["EcO"]),
            "matched_EuE_EcP_pair": str(has_eue and has_ecp),
            "matched_EuE_EcP_AdjP_triplet": str(has_eue and has_ecp and has_adjp),
            "matched_EuE_EcP_EcO_triplet": str(has_eue and has_ecp and has_eco),
        })

    n_scrna = sum(1 for r in sample_rows if r["included_discovery_scrna"] == "True")
    n_pairs = sum(1 for r in donor_rows if r["matched_EuE_EcP_pair"] == "True")
    n_triplets = sum(1 for r in donor_rows if r["matched_EuE_EcP_AdjP_triplet"] == "True")
    assert n_scrna == 31, f"Expected 31 tissue scRNA libraries in GSE179640, got {n_scrna}"
    assert n_pairs == 8, f"Expected 8 matched (EuE, EcP) pairs in GSE179640, got {n_pairs}"
    assert n_triplets == 6, f"Expected 6 matched (EuE, EcP, AdjP) triplets in GSE179640, got {n_triplets}"

    write_tsv(
        ROOT / "manifests/m0_GSE179640_samples.tsv",
        sample_rows,
        list(sample_rows[0].keys()),
    )
    write_tsv(
        ROOT / "manifests/m0_GSE179640_donor_design.tsv",
        donor_rows,
        list(donor_rows[0].keys()),
    )
    return {
        "total_geo_samples": len(sample_rows),
        "tissue_scrna_libraries": n_scrna,
        "donors_total": len(donor_rows),
        "matched_EuE_EcP_pairs": n_pairs,
        "matched_EuE_EcP_AdjP_triplets": n_triplets,
        "ecO_donors": sum(1 for r in donor_rows if r["EcO_gsms"]),
        "ctrl_eue_donors": sum(1 for r in donor_rows if r["Ctrl_EuE_gsms"]),
        "total_scrna_bytes": sum(r["bytes"] for r in sample_rows if r["included_discovery_scrna"] == "True"),
        "skipped_table_lines": skipped,
    }


def audit_gse213216(meta_dir: pathlib.Path, v1_eutopic_ledger: pathlib.Path):
    _, samples, skipped = parse_soft_metadata_only(meta_dir / "GSE213216_family.soft.gz")
    sizes = parse_filelist(meta_dir / "filelist.txt")

    # Read clinical metadata curated in v1's eutopic donor ledger (cycle phase, hormone status, pathology)
    curated_by_gsm = {}
    if v1_eutopic_ledger.exists():
        with v1_eutopic_ledger.open("r", encoding="utf-8", newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["study"] == "GSE213216":
                    curated_by_gsm[r["gsm"]] = r

    all_rows = []
    candidates_by_patient_class = defaultdict(list)

    for s in sorted(samples, key=lambda x: x["gsm"]):
        f = s["fields"]
        c = s["chars"]
        lib_strat = f.get("Sample_library_strategy", [""])[0]
        supp_urls = [u for k, vals in f.items() if k.startswith("Sample_supplementary_file") for u in vals if u != "NONE"]
        supp_names = [u.split("/")[-1] for u in supp_urls]
        bytes_total = sum(sizes.get(fn, 0) for fn in supp_names)
        patient_raw = c.get("patient id", "").strip()
        patient_id = f"P{int(patient_raw):02d}" if patient_raw.isdigit() else "UNLINKED"
        tissue = c.get("tissue", "")
        anat = c.get("anatomic location", "")

        # Deterministic compartment classification from GEO metadata
        site_class = "EXCLUDED"
        excl_reason = ""
        if lib_strat != "RNA-Seq":
            excl_reason = "non_scRNA_spatial_or_other"
        elif not supp_names:
            excl_reason = "no_public_supplementary_matrix_in_GEO"
        elif tissue == "Ovary":
            excl_reason = "whole_ovary_non_endometrial_non_lesion"
        elif tissue == "No endometriosis detected" and "bladder endometriosis" in anat.lower():
            excl_reason = "contradictory_label_no_endo_detected_vs_bladder_endometriosis"
        elif tissue == "Eutopic Endometrium":
            site_class = "EuE"
        elif tissue == "Endometriosis":
            site_class = "EcP"
        elif tissue == "Endometrioma":
            site_class = "EcO"
        elif tissue == "No endometriosis detected":
            site_class = "Unaffected_Peritoneum"
        else:
            excl_reason = f"unrecognized_tissue_{tissue}"

        cur_meta = curated_by_gsm.get(s["gsm"], {})
        row = {
            "gsm": s["gsm"],
            "patient_id": patient_id,
            "patient_num": patient_raw,
            "title": f.get("Sample_title", [""])[0],
            "library_strategy": lib_strat,
            "tissue": tissue,
            "anatomic_location": anat,
            "site_class": site_class,
            "curated_disease": cur_meta.get("disease", "ENDO" if site_class in ("EcP", "EcO", "Unaffected_Peritoneum") else ""),
            "curated_cycle_phase": cur_meta.get("binary_stage", ""),
            "curated_hormonal_treatment": cur_meta.get("hormonal_treatment", ""),
            "supplementary_files": ";".join(supp_names),
            "supplementary_urls": ";".join(supp_urls),
            "bytes": bytes_total,
            "exclusion_reason": excl_reason,
        }
        all_rows.append(row)
        if site_class != "EXCLUDED":
            candidates_by_patient_class[(patient_id, site_class)].append(row)

    # Deterministic non-redundant library selection rule:
    # Select 1 representative library per (patient_id, site_class) by ascending GSM accession.
    # Guarantees 100% coverage of all unique patients in every compartment without wasting
    # bandwidth/compute on redundant multi-biopsy libraries from a single donor.
    selected_gsms = set()
    for key in sorted(candidates_by_patient_class.keys()):
        chosen = sorted(candidates_by_patient_class[key], key=lambda r: r["gsm"])[0]
        selected_gsms.add(chosen["gsm"])

    locked_rows = []
    by_pt_selected = defaultdict(dict)
    for r in all_rows:
        r["selected_for_stage3_download"] = str(r["gsm"] in selected_gsms)
        if r["gsm"] in selected_gsms:
            locked_rows.append(r)
            by_pt_selected[r["patient_id"]][r["site_class"]] = r["gsm"]

    matched_eue_ecp_patients = sorted(
        pt for pt, d in by_pt_selected.items() if "EuE" in d and "EcP" in d
    )
    assert matched_eue_ecp_patients == ["P10", "P11", "P14", "P16"], (
        f"Expected 4 matched (EuE, EcP) patients P10,P11,P14,P16; got {matched_eue_ecp_patients}"
    )

    write_tsv(ROOT / "manifests/m0_GSE213216_samples.tsv", all_rows, list(all_rows[0].keys()))
    write_tsv(ROOT / "manifests/m0_GSE213216_locked_selection.tsv", locked_rows, list(locked_rows[0].keys()))

    comp_counts = Counter(r["site_class"] for r in locked_rows)
    return {
        "total_geo_samples": len(all_rows),
        "locked_libraries_for_download": len(locked_rows),
        "locked_unique_patients": len(by_pt_selected),
        "locked_by_compartment": dict(comp_counts),
        "matched_EuE_EcP_patients": matched_eue_ecp_patients,
        "locked_download_bytes": sum(r["bytes"] for r in locked_rows),
        "skipped_table_lines": skipped,
    }


def audit_gse141549(meta_dir: pathlib.Path):
    _, samples, skipped = parse_soft_metadata_only(meta_dir / "GSE141549_family.soft.gz")
    wb = openpyxl.load_workbook(meta_dir / "GSE141549_Sample_link.xlsx", read_only=True, data_only=True)
    ws = wb["Sheet1"]
    title_to_subject = {}
    for r in list(ws.iter_rows(values_only=True))[1:]:
        subj = str(r[0]).strip() if r[0] else ""
        if not subj:
            continue
        for cell in r[1:]:
            if cell:
                t = str(cell).strip()
                title_to_subject[t] = subj

    sample_rows = []
    by_subject = defaultdict(lambda: defaultdict(list))

    for s in sorted(samples, key=lambda x: x["gsm"]):
        f = s["fields"]
        c = s["chars"]
        title = f.get("Sample_title", [""])[0].strip()
        desc = f.get("Sample_description", [""])[0].strip()
        platform = f.get("Sample_platform_id", [""])[0].strip()
        subject_id = title_to_subject.get(title, "")
        if not subject_id:
            raise ValueError(f"Unlinked sample in GSE141549: {s['gsm']} {title}")

        parts = title.split()
        raw_code = parts[2] if len(parts) >= 3 else title
        is_tech_rep = "replicate" in title.lower() or "replicate" in c.get("tissue", "").lower()

        if raw_code in ("PeLR", "PeLB", "PeLW"):
            major_group = "Peritoneal_Lesion"
        elif raw_code in ("DiEIn", "DiEB", "REV", "SuL"):
            major_group = "DIE"
        elif raw_code == "OMA":
            major_group = "OMA"
        elif raw_code in ("PE", "PP", "CE", "CP"):
            major_group = raw_code
        else:
            raise ValueError(f"Unexpected GSE141549 tissue code: {raw_code}")

        age_raw = c.get("age", "").strip()
        stage_raw = c.get("disease stage", "").strip()
        phase_raw = c.get("cycle phase", "unknown").strip().lower()
        if phase_raw in ("proliferative", "secretory", "menstruation", "medication"):
            phase_clean = phase_raw
        else:
            phase_clean = "unknown"

        disease_status = "Control" if (raw_code in ("CE", "CP") or stage_raw.lower() == "healthy") else "Endometriosis"

        sample_rows.append({
            "gsm": s["gsm"],
            "subject_id": subject_id,
            "sample_title": title,
            "platform_id": platform,
            "tissue_code": raw_code,
            "major_group": major_group,
            "disease_status": disease_status,
            "disease_stage": stage_raw or "NA",
            "cycle_phase": phase_clean,
            "age": age_raw or "NA",
            "is_technical_replicate": str(is_tech_rep),
            "description": desc,
        })
        by_subject[subject_id][major_group].append(s["gsm"])

    subject_rows = []
    for subj in sorted(by_subject.keys()):
        g = by_subject[subj]
        subject_rows.append({
            "subject_id": subj,
            "CE_n": len(g["CE"]),
            "CP_n": len(g["CP"]),
            "PE_n": len(g["PE"]),
            "PP_n": len(g["PP"]),
            "Peritoneal_Lesion_n": len(g["Peritoneal_Lesion"]),
            "DIE_n": len(g["DIE"]),
            "OMA_n": len(g["OMA"]),
            "has_paired_PE_Peritoneal": str(len(g["PE"]) > 0 and len(g["Peritoneal_Lesion"]) > 0),
            "has_paired_PP_Peritoneal": str(len(g["PP"]) > 0 and len(g["Peritoneal_Lesion"]) > 0),
            "has_triplet_PE_PP_Peritoneal": str(len(g["PE"]) > 0 and len(g["PP"]) > 0 and len(g["Peritoneal_Lesion"]) > 0),
            "has_paired_PP_DIE": str(len(g["PP"]) > 0 and len(g["DIE"]) > 0),
        })

    assert len(sample_rows) == 408, f"Expected 408 samples in GSE141549, got {len(sample_rows)}"
    assert len(subject_rows) == 162, f"Expected 162 linked subjects in GSE141549, got {len(subject_rows)}"

    write_tsv(ROOT / "manifests/m0_GSE141549_samples.tsv", sample_rows, list(sample_rows[0].keys()))
    write_tsv(ROOT / "manifests/m0_GSE141549_subjects.tsv", subject_rows, list(subject_rows[0].keys()))

    return {
        "total_samples": len(sample_rows),
        "unique_subjects": len(subject_rows),
        "samples_by_major_group": dict(Counter(r["major_group"] for r in sample_rows)),
        "subjects_by_major_group": {
            mg: sum(1 for subj in by_subject.values() if len(subj[mg]) > 0)
            for mg in ["CE", "CP", "PE", "PP", "Peritoneal_Lesion", "DIE", "OMA"]
        },
        "paired_PE_Peritoneal_subjects": sum(1 for r in subject_rows if r["has_paired_PE_Peritoneal"] == "True"),
        "paired_PP_Peritoneal_subjects": sum(1 for r in subject_rows if r["has_paired_PP_Peritoneal"] == "True"),
        "triplet_PE_PP_Peritoneal_subjects": sum(1 for r in subject_rows if r["has_triplet_PE_PP_Peritoneal"] == "True"),
        "paired_PP_DIE_subjects": sum(1 for r in subject_rows if r["has_paired_PP_DIE"] == "True"),
        "skipped_table_lines": skipped,
    }


def build_genetics_prior(gen_dir: pathlib.Path, hgnc: dict):
    # 1. Parse Open Targets MONDO_0005133 (EFO:0001065) target pages
    ot_by_ensg = {}
    meta_info = None
    for p in sorted(gen_dir.glob("ot_targets_page_*.json")):
        d = json.loads(p.read_text(encoding="utf-8"))
        if meta_info is None:
            meta_info = d.get("meta")
        for r in d["disease"]["associatedTargets"]["rows"]:
            t = r["target"]
            eg = t["id"]
            if not ENSG_RE.fullmatch(eg):
                continue
            dt_scores = {x["id"]: float(x["score"]) for x in r.get("datatypeScores", [])}
            ds_scores = {x["id"]: float(x["score"]) for x in r.get("datasourceScores", [])}
            tract = t.get("tractability", []) or []
            sm_clin = any(x.get("modality") == "SM" and x.get("value") and "Clinical" in x.get("label", "") or (x.get("modality") == "SM" and x.get("value") and x.get("label") == "Approved Drug") for x in tract)
            sm_pocket = any(x.get("modality") == "SM" and x.get("value") for x in tract)
            ab_clin = any(x.get("modality") == "AB" and x.get("value") and ("Clinical" in x.get("label", "") or x.get("label") == "Approved Drug") for x in tract)
            ab_tract = any(x.get("modality") == "AB" and x.get("value") for x in tract)
            protac_tract = any(x.get("modality") == "PR" and x.get("value") for x in tract)

            prev = ot_by_ensg.get(eg)
            if prev is None or float(r["score"]) > prev["ot_overall_score"]:
                ot_by_ensg[eg] = {
                    "ensembl_gene_id": eg,
                    "symbol": hgnc["ensg_to_symbol"].get(eg, t.get("approvedSymbol", "")),
                    "ot_overall_score": float(r["score"]),
                    "ot_genetic_association": dt_scores.get("genetic_association", 0.0),
                    "ot_l2g_credible_set": ds_scores.get("gwas_credible_sets", 0.0),
                    "ot_somatic_mutation": dt_scores.get("somatic_mutation", 0.0),
                    "ot_clinical_drug": dt_scores.get("clinical", 0.0),
                    "ot_rna_expression": dt_scores.get("rna_expression", 0.0),
                    "ot_animal_model": dt_scores.get("animal_model", 0.0),
                    "tractable_sm_clinical": sm_clin,
                    "tractable_sm_any": sm_pocket,
                    "tractable_ab_clinical": ab_clin,
                    "tractable_ab_any": ab_tract,
                    "tractable_protac": protac_tract,
                }

    # 2. Parse Rahmioglu et al. 2023 Nat Genet supplement (41588_2023_1323_MOESM4_ESM.xlsx)
    wb = openpyxl.load_workbook(gen_dir / "41588_2023_1323_MOESM4_ESM.xlsx", read_only=True, data_only=True)

    def get_sheet(prefix: str):
        return wb[next(n for n in wb.sheetnames if n.startswith(prefix + "."))]

    def extract_symbols(ws, header_row: int, col_name: str) -> set[str]:
        out = set()
        hdr = None
        for idx_r, r in enumerate(ws.iter_rows(values_only=True), 1):
            if idx_r < header_row:
                continue
            if idx_r == header_row:
                hdr = [str(x).strip() if x is not None else "" for x in r]
                col_i = hdr.index(col_name)
                continue
            if col_i >= len(r) or r[col_i] is None:
                continue
            val = str(r[col_i]).replace("\xa0", " ").strip()
            for tok in re.split(r"[,;/]", val):
                tok = re.sub(r"[*^]+$", "", tok.strip())
                if tok and tok != "." and re.fullmatch(r"[A-Za-z][A-Za-z0-9._-]{0,24}", tok):
                    out.add(tok)
        return out

    # Supp3 Locus column (42 GWS lead loci)
    ws_s3 = get_sheet("Supp3")
    gws_loci_syms = set()
    for idx_r, r in enumerate(ws_s3.iter_rows(values_only=True), 1):
        if idx_r <= 6 or not r or r[0] is None or r[2] is None:
            continue
        loc_str = str(r[0]).replace("\xa0", " ").strip()
        if not str(r[2]).strip().startswith("rs"):
            continue
        for tok in re.split(r"[/,; ]+", loc_str):
            tok = tok.strip()
            if tok and re.fullmatch(r"[A-Za-z][A-Za-z0-9._-]{1,24}", tok):
                gws_loci_syms.add(tok)

    cred_nearest_syms = extract_symbols(get_sheet("Supp11"), 2, "Nearest gene")
    coding_cred_syms = extract_symbols(get_sheet("Supp12"), 2, "Gene")
    smr_uterus_syms = extract_symbols(get_sheet("Supp14"), 2, "Gene")
    smr_endo_exp_syms = extract_symbols(get_sheet("Supp15"), 2, "Gene")
    smr_blood_exp_syms = extract_symbols(get_sheet("Supp16"), 2, "Gene")
    smr_endo_meth_syms = extract_symbols(get_sheet("Supp17"), 2, "Closest Gene")
    smr_blood_meth_syms = extract_symbols(get_sheet("Supp18"), 2, "Closest Gene")

    rah_sets = {
        "rahmioglu_gws42": set(),
        "rahmioglu_cred_nearest": set(),
        "rahmioglu_coding_cred": set(),
        "rahmioglu_smr_uterus_endometrium": set(),
        "rahmioglu_smr_any": set(),
    }
    for sym in gws_loci_syms:
        eg, _ = resolve_symbol(sym, hgnc)
        if eg:
            rah_sets["rahmioglu_gws42"].add(eg)
    for sym in cred_nearest_syms:
        eg, _ = resolve_symbol(sym, hgnc)
        if eg:
            rah_sets["rahmioglu_cred_nearest"].add(eg)
    for sym in coding_cred_syms:
        eg, _ = resolve_symbol(sym, hgnc)
        if eg:
            rah_sets["rahmioglu_coding_cred"].add(eg)
    for sym in (smr_uterus_syms | smr_endo_exp_syms | smr_endo_meth_syms):
        eg, _ = resolve_symbol(sym, hgnc)
        if eg:
            rah_sets["rahmioglu_smr_uterus_endometrium"].add(eg)
    for sym in (smr_uterus_syms | smr_endo_exp_syms | smr_blood_exp_syms | smr_endo_meth_syms | smr_blood_meth_syms):
        eg, _ = resolve_symbol(sym, hgnc)
        if eg:
            rah_sets["rahmioglu_smr_any"].add(eg)

    all_prior_genes = sorted(
        set(ot_by_ensg.keys())
        | rah_sets["rahmioglu_gws42"]
        | rah_sets["rahmioglu_cred_nearest"]
        | rah_sets["rahmioglu_coding_cred"]
        | rah_sets["rahmioglu_smr_any"]
    )

    out_rows = []
    for eg in all_prior_genes:
        ot = ot_by_ensg.get(eg, {})
        ot_ga = float(ot.get("ot_genetic_association", 0.0))
        in_gwas = eg in rah_sets["rahmioglu_cred_nearest"] or eg in rah_sets["rahmioglu_gws42"]
        in_coding = eg in rah_sets["rahmioglu_coding_cred"]
        in_smr_ut = eg in rah_sets["rahmioglu_smr_uterus_endometrium"]
        in_smr_any = eg in rah_sets["rahmioglu_smr_any"]

        # Pre-specified bounded independent human genetic causality prior in [0, 1]
        # Uses ONLY human germline genetics / fine-mapped L2G / SMR-eQTL / coding credible variants
        # (never literature text-mining or cohort RNA expression!)
        g_prior = min(
            1.0,
            0.60 * ot_ga
            + 0.35 * (1.0 if in_gwas else 0.0)
            + 0.45 * (1.0 if in_smr_any else 0.0)
            + 0.25 * (1.0 if in_coding else 0.0),
        )

        out_rows.append({
            "ensembl_gene_id": eg,
            "symbol": ot.get("symbol") or hgnc["ensg_to_symbol"].get(eg, ""),
            "genetic_prior_score": round(g_prior, 6),
            "has_human_genetic_support": str(g_prior > 0.0),
            "ot_overall_score": round(float(ot.get("ot_overall_score", 0.0)), 6),
            "ot_genetic_association": round(ot_ga, 6),
            "ot_l2g_credible_set": round(float(ot.get("ot_l2g_credible_set", 0.0)), 6),
            "ot_somatic_mutation": round(float(ot.get("ot_somatic_mutation", 0.0)), 6),
            "ot_clinical_drug": round(float(ot.get("ot_clinical_drug", 0.0)), 6),
            "rahmioglu_gws42": str(eg in rah_sets["rahmioglu_gws42"]),
            "rahmioglu_cred_nearest": str(eg in rah_sets["rahmioglu_cred_nearest"]),
            "rahmioglu_coding_cred": str(in_coding),
            "rahmioglu_smr_uterus_endo": str(in_smr_ut),
            "rahmioglu_smr_any": str(in_smr_any),
            "tractable_sm_clinical": str(bool(ot.get("tractable_sm_clinical", False))),
            "tractable_sm_any": str(bool(ot.get("tractable_sm_any", False))),
            "tractable_ab_clinical": str(bool(ot.get("tractable_ab_clinical", False))),
            "tractable_ab_any": str(bool(ot.get("tractable_ab_any", False))),
            "tractable_protac": str(bool(ot.get("tractable_protac", False))),
        })

    out_rows.sort(key=lambda r: (-r["genetic_prior_score"], -r["ot_overall_score"], r["ensembl_gene_id"]))
    write_tsv(ROOT / "manifests/human_genetics_prior_v2.tsv", out_rows, list(out_rows[0].keys()))

    return {
        "ot_meta": meta_info,
        "ot_unique_targets": len(ot_by_ensg),
        "ot_genetic_targets_gt0": sum(1 for r in out_rows if r["ot_genetic_association"] > 0),
        "ot_genetic_targets_ge0_5": sum(1 for r in out_rows if r["ot_genetic_association"] >= 0.5),
        "rahmioglu_gws42_genes": len(rah_sets["rahmioglu_gws42"]),
        "rahmioglu_cred_nearest_genes": len(rah_sets["rahmioglu_cred_nearest"]),
        "rahmioglu_coding_cred_genes": len(rah_sets["rahmioglu_coding_cred"]),
        "rahmioglu_smr_uterus_endo_genes": len(rah_sets["rahmioglu_smr_uterus_endometrium"]),
        "rahmioglu_smr_any_genes": len(rah_sets["rahmioglu_smr_any"]),
        "combined_genetic_support_genes": sum(1 for r in out_rows if r["genetic_prior_score"] > 0),
    }


def build_canonical_markers_and_lr(hgnc: dict):
    # De novo canonical lineage marker panel (7 compartments; zero HECA transferred labels)
    panel_spec = [
        ("stromal_fibroblast", ["DCN", "LUM", "PDGFRA", "COL1A1", "COL1A2", "COL3A1", "COL6A1", "MMP2"],
         "Canonical pan-stromal ECM/PDGFRA+ fibroblast identity"),
        ("perivascular_mural", ["RGS5", "NOTCH3", "MYH11", "MCAM", "TAGLN", "ACTA2", "CNN1", "CSPG4"],
         "Canonical pericyte and vascular smooth muscle mural identity"),
        ("epithelial", ["EPCAM", "KRT8", "KRT18", "PAX8", "WFDC2", "MUC1", "CDH1", "FOXJ1"],
         "Canonical Mullerian glandular/luminal/ciliated epithelial identity"),
        ("endothelial", ["PECAM1", "VWF", "CLDN5", "CDH5", "KDR", "ACKR1", "ERG", "EMCN"],
         "Canonical vascular endothelial identity"),
        ("myeloid_macrophage", ["LYZ", "AIF1", "CD68", "C1QA", "C1QB", "TYROBP", "FCER1G", "CD14"],
         "Canonical monocyte/macrophage/dendritic myeloid identity"),
        ("lymphoid", ["CD3D", "CD3E", "TRBC1", "IL7R", "NKG7", "GNLY", "MS4A1", "CD79A"],
         "Canonical T/NK/B lymphoid identity"),
        ("mesothelial", ["MSLN", "UPK3B", "CALB2", "LRRN4", "ITLN1", "PRG4", "KRT5", "BNC1"],
         "Canonical serosal peritoneal mesothelial identity (ambient host-tissue audit & isolation)"),
    ]
    marker_rows = []
    for lineage, syms, rationale in panel_spec:
        for sym in syms:
            eg, how = resolve_symbol(sym, hgnc)
            if not eg:
                raise ValueError(f"Canonical marker {sym} failed HGNC resolution")
            marker_rows.append({
                "lineage": lineage,
                "symbol": sym,
                "ensembl_gene_id": eg,
                "mapping_method": how,
                "rationale": rationale,
            })
    write_tsv(ROOT / "manifests/canonical_lineage_markers.tsv", marker_rows, list(marker_rows[0].keys()))

    # Curated canonical Ligand-Receptor interactions across vascular/myeloid <-> stromal niches
    # Sourced from canonical human L-R literature (Ramilowski 2015 / CellPhoneDB / NicheNet core axes)
    lr_raw = [
        # Angiogenic & endothelial-stromal crosstalk
        ("VEGFA", "KDR", "VEGF_angiogenesis"), ("VEGFA", "FLT1", "VEGF_angiogenesis"),
        ("VEGFA", "NRP1", "VEGF_angiogenesis"), ("VEGFA", "NRP2", "VEGF_angiogenesis"),
        ("VEGFB", "FLT1", "VEGF_angiogenesis"), ("VEGFC", "FLT4", "VEGF_angiogenesis"),
        ("VEGFC", "KDR", "VEGF_angiogenesis"), ("PGF", "FLT1", "VEGF_angiogenesis"),
        ("ANGPT1", "TEK", "Angiopoietin_Tie"), ("ANGPT2", "TEK", "Angiopoietin_Tie"),
        ("ANGPT2", "ITGB1", "Angiopoietin_Tie"), ("PDGFA", "PDGFRA", "PDGF_stromal_vascular"),
        ("PDGFB", "PDGFRB", "PDGF_stromal_vascular"), ("PDGFB", "PDGFRA", "PDGF_stromal_vascular"),
        ("PDGFC", "PDGFRA", "PDGF_stromal_vascular"), ("PDGFD", "PDGFRB", "PDGF_stromal_vascular"),
        ("JAG1", "NOTCH1", "Notch_vascular_niche"), ("JAG1", "NOTCH2", "Notch_vascular_niche"),
        ("JAG1", "NOTCH3", "Notch_vascular_niche"), ("JAG2", "NOTCH1", "Notch_vascular_niche"),
        ("DLL1", "NOTCH1", "Notch_vascular_niche"), ("DLL4", "NOTCH1", "Notch_vascular_niche"),
        ("DLL4", "NOTCH4", "Notch_vascular_niche"), ("EFNB2", "EPHB4", "Ephrin_vascular"),
        ("EFNB1", "EPHB2", "Ephrin_vascular"), ("EFNA1", "EPHA2", "Ephrin_vascular"),
        ("EDN1", "EDNRA", "Endothelin_vascular"), ("EDN1", "EDNRB", "Endothelin_vascular"),
        ("ADM", "CALCRL", "Adrenomedullin_vascular"), ("ADM", "RAMP2", "Adrenomedullin_vascular"),
        ("SLIT2", "ROBO1", "Slit_Robo_guidance"), ("SLIT2", "ROBO4", "Slit_Robo_guidance"),
        ("SLIT3", "ROBO4", "Slit_Robo_guidance"), ("SEMA3A", "NRP1", "Semaphorin_guidance"),
        ("SEMA3C", "NRP2", "Semaphorin_guidance"), ("SEMA3C", "PLXND1", "Semaphorin_guidance"),
        ("SEMA3F", "NRP2", "Semaphorin_guidance"), ("ACKR1", "CXCL8", "DARC_endothelial_chemokine"),
        ("ACKR1", "CCL2", "DARC_endothelial_chemokine"), ("ACKR1", "CXCL1", "DARC_endothelial_chemokine"),
        ("ACKR1", "CXCL2", "DARC_endothelial_chemokine"), ("ACKR1", "CXCL6", "DARC_endothelial_chemokine"),
        # Myeloid-stromal inflammatory cytokines & growth factors
        ("IL1A", "IL1R1", "IL1_inflammatory_niche"), ("IL1A", "IL1RAP", "IL1_inflammatory_niche"),
        ("IL1B", "IL1R1", "IL1_inflammatory_niche"), ("IL1B", "IL1RAP", "IL1_inflammatory_niche"),
        ("IL6", "IL6R", "IL6_gp130_jak_stat"), ("IL6", "IL6ST", "IL6_gp130_jak_stat"),
        ("OSM", "OSMR", "IL6_gp130_jak_stat"), ("OSM", "IL6ST", "IL6_gp130_jak_stat"),
        ("LIF", "LIFR", "IL6_gp130_jak_stat"), ("LIF", "IL6ST", "IL6_gp130_jak_stat"),
        ("TNF", "TNFRSF1A", "TNF_inflammatory"), ("TNF", "TNFRSF1B", "TNF_inflammatory"),
        ("TNFSF12", "TNFRSF12A", "TWEAK_Fn14_remodeling"), ("TNFSF13B", "TNFRSF13C", "BAFF_immune"),
        ("CSF1", "CSF1R", "Macrophage_trophic_CSF1"), ("IL34", "CSF1R", "Macrophage_trophic_CSF1"),
        ("GAS6", "AXL", "TAM_effero_fibrosis"), ("GAS6", "MERTK", "TAM_effero_fibrosis"),
        ("GAS6", "TYRO3", "TAM_effero_fibrosis"), ("PROS1", "AXL", "TAM_effero_fibrosis"),
        ("PROS1", "MERTK", "TAM_effero_fibrosis"), ("MIF", "CD74", "MIF_macrophage"),
        ("MIF", "CXCR4", "MIF_macrophage"), ("C3", "C3AR1", "Complement_myeloid"),
        ("C3", "ITGAM", "Complement_myeloid"), ("C3", "ITGB2", "Complement_myeloid"),
        ("C5", "C5AR1", "Complement_myeloid"), ("RARRES2", "CMKLR1", "Chemerin_myeloid"),
        # TGF-beta / BMP / Wnt / Morphogen fibrotic remodeling
        ("TGFB1", "TGFBR1", "TGFB_fibrosis"), ("TGFB1", "TGFBR2", "TGFB_fibrosis"),
        ("TGFB1", "TGFBR3", "TGFB_fibrosis"), ("TGFB2", "TGFBR2", "TGFB_fibrosis"),
        ("TGFB3", "TGFBR2", "TGFB_fibrosis"), ("INHBA", "ACVR2A", "Activin_TGFB"),
        ("INHBA", "ACVR1B", "Activin_TGFB"), ("BMP2", "BMPR1A", "BMP_morphogen"),
        ("BMP2", "BMPR2", "BMP_morphogen"), ("BMP4", "BMPR1A", "BMP_morphogen"),
        ("BMP4", "BMPR2", "BMP_morphogen"), ("BMP6", "ACVR1", "BMP_morphogen"),
        ("GREM1", "BMPR2", "BMP_antagonist_niche"), ("WNT4", "FZD2", "Wnt_Mullerian_remodeling"),
        ("WNT4", "FZD6", "Wnt_Mullerian_remodeling"), ("WNT4", "LRP6", "Wnt_Mullerian_remodeling"),
        ("WNT5A", "FZD5", "Wnt_Mullerian_remodeling"), ("WNT5A", "ROR2", "Wnt_Mullerian_remodeling"),
        ("WNT2", "FZD4", "Wnt_Mullerian_remodeling"), ("RSPO3", "LGR4", "RSPO_Wnt_amplifier"),
        ("RSPO3", "LGR5", "RSPO_Wnt_amplifier"), ("RSPO1", "LGR4", "RSPO_Wnt_amplifier"),
        ("SFRP2", "FZD2", "Wnt_modulator"),
        # Chemokine leukocyte recruitment & stromal retention
        ("CCL2", "CCR2", "Chemokine_recruitment"), ("CCL3", "CCR1", "Chemokine_recruitment"),
        ("CCL3", "CCR5", "Chemokine_recruitment"), ("CCL4", "CCR5", "Chemokine_recruitment"),
        ("CCL5", "CCR5", "Chemokine_recruitment"), ("CCL5", "CCR1", "Chemokine_recruitment"),
        ("CCL7", "CCR2", "Chemokine_recruitment"), ("CCL8", "CCR2", "Chemokine_recruitment"),
        ("CCL14", "CCR1", "Chemokine_recruitment"), ("CCL19", "CCR7", "Chemokine_recruitment"),
        ("CCL21", "CCR7", "Chemokine_recruitment"), ("CXCL12", "CXCR4", "CXCL12_stromal_niche"),
        ("CXCL12", "ACKR3", "CXCL12_stromal_niche"), ("CXCL14", "CXCR4", "CXCL12_stromal_niche"),
        ("CXCL8", "CXCR1", "Chemokine_recruitment"), ("CXCL8", "CXCR2", "Chemokine_recruitment"),
        ("CXCL1", "CXCR2", "Chemokine_recruitment"), ("CXCL2", "CXCR2", "Chemokine_recruitment"),
        ("CXCL9", "CXCR3", "Chemokine_recruitment"), ("CXCL10", "CXCR3", "Chemokine_recruitment"),
        ("CXCL16", "CXCR6", "Chemokine_recruitment"),
        # RTK trophic & neuro-angiogenic pain/growth signaling
        ("IGF1", "IGF1R", "IGF_FGF_HGF_trophic"), ("IGF2", "IGF1R", "IGF_FGF_HGF_trophic"),
        ("IGF2", "IGF2R", "IGF_FGF_HGF_trophic"), ("HGF", "MET", "IGF_FGF_HGF_trophic"),
        ("FGF2", "FGFR1", "IGF_FGF_HGF_trophic"), ("FGF7", "FGFR2", "IGF_FGF_HGF_trophic"),
        ("FGF9", "FGFR1", "IGF_FGF_HGF_trophic"), ("FGF9", "FGFR3", "IGF_FGF_HGF_trophic"),
        ("HBEGF", "EGFR", "EGF_family"), ("AREG", "EGFR", "EGF_family"),
        ("EREG", "EGFR", "EGF_family"), ("NGF", "NTRK1", "Neurotrophin_pain_niche"),
        ("NGF", "NGFR", "Neurotrophin_pain_niche"), ("BDNF", "NTRK2", "Neurotrophin_pain_niche"),
        # ECM-Integrin / CD44 adhesion & mechanotransduction
        ("FN1", "ITGA5", "ECM_integrin_adhesion"), ("FN1", "ITGB1", "ECM_integrin_adhesion"),
        ("FN1", "ITGAV", "ECM_integrin_adhesion"), ("FN1", "ITGB3", "ECM_integrin_adhesion"),
        ("FN1", "CD44", "ECM_integrin_adhesion"), ("FN1", "SDC4", "ECM_integrin_adhesion"),
        ("SPP1", "CD44", "SPP1_macrophage_adhesion"), ("SPP1", "ITGAV", "SPP1_macrophage_adhesion"),
        ("SPP1", "ITGB1", "SPP1_macrophage_adhesion"), ("SPP1", "ITGB5", "SPP1_macrophage_adhesion"),
        ("COL1A1", "ITGA1", "ECM_integrin_adhesion"), ("COL1A1", "ITGA2", "ECM_integrin_adhesion"),
        ("COL1A1", "ITGB1", "ECM_integrin_adhesion"), ("COL1A2", "ITGB1", "ECM_integrin_adhesion"),
        ("COL3A1", "ITGA1", "ECM_integrin_adhesion"), ("COL3A1", "ITGB1", "ECM_integrin_adhesion"),
        ("COL4A1", "ITGA1", "ECM_integrin_adhesion"), ("COL4A1", "ITGB1", "ECM_integrin_adhesion"),
        ("COL6A1", "ITGB1", "ECM_integrin_adhesion"), ("COL6A2", "CD44", "ECM_integrin_adhesion"),
        ("THBS1", "CD36", "Thrombospondin_adhesion"), ("THBS1", "CD47", "Thrombospondin_adhesion"),
        ("THBS1", "ITGAV", "Thrombospondin_adhesion"), ("THBS2", "CD36", "Thrombospondin_adhesion"),
        ("TNC", "ITGAV", "ECM_integrin_adhesion"), ("TNC", "ITGB1", "ECM_integrin_adhesion"),
        ("LAMA4", "ITGA6", "ECM_integrin_adhesion"), ("LAMB1", "ITGB1", "ECM_integrin_adhesion"),
        ("CCN2", "ITGAV", "Matricellular_CCN"), ("CCN2", "LRP1", "Matricellular_CCN"),
        ("CCN1", "ITGAV", "Matricellular_CCN"), ("LGALS1", "PTPRC", "Galectin_immune"),
        ("LGALS3", "LRP1", "Galectin_immune"), ("LGALS9", "HAVCR2", "Galectin_immune"),
        ("ICAM1", "ITGAL", "Endothelial_leukocyte_adhesion"), ("ICAM1", "ITGAM", "Endothelial_leukocyte_adhesion"),
        ("VCAM1", "ITGA4", "Endothelial_leukocyte_adhesion"), ("VCAM1", "ITGB1", "Endothelial_leukocyte_adhesion"),
        ("SELE", "SELPLG", "Endothelial_leukocyte_adhesion"), ("SELP", "SELPLG", "Endothelial_leukocyte_adhesion"),
    ]

    lr_rows = []
    seen_pairs = set()
    for lig, rec, pathway in lr_raw:
        eg_l, _ = resolve_symbol(lig, hgnc)
        eg_r, _ = resolve_symbol(rec, hgnc)
        if not eg_l or not eg_r or eg_l == eg_r:
            raise ValueError(f"Curated L-R pair ({lig}, {rec}) failed unique HGNC->ENSG resolution")
        if (eg_l, eg_r) in seen_pairs:
            continue
        seen_pairs.add((eg_l, eg_r))
        lr_rows.append({
            "pair_id": f"{lig}_{rec}",
            "ligand_symbol": lig,
            "ligand_ensembl_id": eg_l,
            "receptor_symbol": rec,
            "receptor_ensembl_id": eg_r,
            "pathway_group": pathway,
        })
    write_tsv(ROOT / "manifests/curated_niche_lr_pairs.tsv", lr_rows, list(lr_rows[0].keys()))
    return {
        "canonical_markers_count": len(marker_rows),
        "lineages_count": len(panel_spec),
        "curated_lr_pairs_count": len(lr_rows),
        "curated_lr_unique_genes": len({r["ligand_ensembl_id"] for r in lr_rows} | {r["receptor_ensembl_id"] for r in lr_rows}),
    }


def build_string_gene_network(ref_dir: pathlib.Path, hgnc: dict):
    """Deterministically build STRING v12.5 direct (experiments/database >= 700) Ensembl gene graph."""
    info_path = ref_dir / "9606.protein.info.v12.5.txt.gz"
    links_path = ref_dir / "9606.protein.links.full.v12.5.txt.gz"

    # Map protein_id -> unique stable ENSG following strict approved/unique-previous symbol rule
    prev_multi = defaultdict(list)
    for sym, recs in hgnc["approved_multi"].items():
        for r in recs:
            for p in filter(None, r.get("prev_symbol", "").split("|")):
                prev_multi[p.strip()].append(r)

    prot_to_ensg = {}
    with gzip.open(info_path, "rt", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            pid = r["#string_protein_id"]
            name = r["preferred_name"]
            recs = hgnc["approved_multi"].get(name) or prev_multi.get(name, [])
            hgnc_ids = {x["hgnc_id"] for x in recs}
            targets = {x["ensembl_gene_id"] for x in recs if ENSG_RE.fullmatch(x["ensembl_gene_id"])}
            if len(hgnc_ids) == 1 and len(targets) == 1:
                prot_to_ensg[pid] = next(iter(targets))

    gene_edges = {}
    source_rows = 0
    qualifying_protein_rows = 0
    with gzip.open(links_path, "rt", encoding="utf-8") as f:
        header = next(f).split()
        i_exp = header.index("experiments")
        i_db = header.index("database")
        for line in f:
            source_rows += 1
            parts = line.split()
            score = max(int(parts[i_exp]), int(parts[i_db]))
            if score < 700:
                continue
            qualifying_protein_rows += 1
            ga = prot_to_ensg.get(parts[0])
            gb = prot_to_ensg.get(parts[1])
            if not ga or not gb or ga == gb:
                continue
            key = (ga, gb) if ga < gb else (gb, ga)
            if score > gene_edges.get(key, 0):
                gene_edges[key] = score

    out_path = ROOT / "data/interim/string_v12_5_direct_gene_edges.tsv.gz"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("wb") as raw_f, gzip.GzipFile(filename="", mode="wb", fileobj=raw_f, mtime=0) as gz:
        gz.write(b"gene1\tgene2\tweight\n")
        for (ga, gb), sc in sorted(gene_edges.items()):
            gz.write(f"{ga}\t{gb}\t{sc / 1000.0:.3f}\n".encode("utf-8"))

    nodes = {g for pair in gene_edges for g in pair}
    summary = {
        "string_version": "12.5",
        "species_taxid": 9606,
        "channels_permitted": ["experiments", "database"],
        "minimum_score": 700,
        "source_protein_rows": source_rows,
        "qualifying_directed_protein_rows": qualifying_protein_rows,
        "mapped_proteins": len(prot_to_ensg),
        "undirected_gene_edges": len(gene_edges),
        "unique_gene_nodes": len(nodes),
        "output_relpath": "data/interim/string_v12_5_direct_gene_edges.tsv.gz",
        "output_sha256": sha256_file(out_path),
        "output_bytes": out_path.stat().st_size,
    }
    (ROOT / "manifests/string_network_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main():
    raw = ROOT / "data/raw"
    hgnc_path = raw / "reference/hgnc_complete_set.tsv"
    hgnc = load_hgnc(hgnc_path)

    v1_ledger = pathlib.Path(
        "/home/seedna2/.gemini/antigravity/brain/c94c7cdc-dbd8-40be-8165-c1c697fb97c6"
        "/scratch/repos/endotarget-evidence/manifests/cohort_extension_eutopic_donor_ledger.tsv"
    )

    g179 = audit_gse179640(raw / "metadata/GSE179640")
    g213 = audit_gse213216(raw / "metadata/GSE213216", v1_ledger)
    g141 = audit_gse141549(raw / "metadata/GSE141549")
    gen = build_genetics_prior(raw / "genetics", hgnc)
    mlr = build_canonical_markers_and_lr(hgnc)
    snet = build_string_gene_network(raw / "reference", hgnc)

    tracked_raw_files = [
        raw / "metadata/GSE179640/GSE179640_family.soft.gz",
        raw / "metadata/GSE179640/filelist.txt",
        raw / "metadata/GSE213216/GSE213216_family.soft.gz",
        raw / "metadata/GSE213216/filelist.txt",
        raw / "metadata/GSE141549/GSE141549_family.soft.gz",
        raw / "metadata/GSE141549/filelist.txt",
        raw / "metadata/GSE141549/GSE141549_Sample_link.xlsx",
        raw / "genetics/41588_2023_1323_MOESM4_ESM.xlsx",
        raw / "genetics/ot_gwas_studies.json",
        raw / "genetics/ot_disease_drugs.json",
        *[raw / f"genetics/ot_targets_page_{i}.json" for i in range(8)],
        raw / "reference/hgnc_complete_set.tsv",
        raw / "reference/9606.protein.info.v12.5.txt.gz",
        raw / "reference/9606.protein.links.full.v12.5.txt.gz",
    ]
    committed_manifests = sorted((ROOT / "manifests").glob("*"))

    receipt = {
        "schema_version": 1,
        "stage": "M0_METADATA_AND_GENETICS_AUDIT",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "attestations": {
            "cohort_expression_matrices_downloaded": False,
            "cohort_expression_values_read": False,
            "heca_cell_annotations_used": False,
            "gse141549_soft_tables_skipped": True,
        },
        "cohorts": {
            "GSE179640_discovery": g179,
            "GSE213216_external_scrna": g213,
            "GSE141549_external_bulk": g141,
        },
        "genetics_and_priors": gen,
        "markers_and_lr": mlr,
        "string_network": snet,
        "raw_input_sha256": {
            str(p.relative_to(ROOT)): {"sha256": sha256_file(p), "bytes": p.stat().st_size}
            for p in tracked_raw_files
        },
        "manifest_sha256": {
            str(p.relative_to(ROOT)): {"sha256": sha256_file(p), "bytes": p.stat().st_size}
            for p in committed_manifests
        },
    }
    rec_path = ROOT / "receipts/m0_inputs_and_priors.json"
    rec_path.parent.mkdir(parents=True, exist_ok=True)
    rec_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "GSE179640": g179,
        "GSE213216": g213,
        "GSE141549": g141,
        "genetics": gen,
        "markers_lr": mlr,
        "string": {k: snet[k] for k in ("undirected_gene_edges", "unique_gene_nodes", "output_sha256")},
    }, indent=2))


if __name__ == "__main__":
    main()
