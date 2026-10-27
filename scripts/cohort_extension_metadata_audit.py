#!/usr/bin/env python3
"""Metadata-only feasibility audit of two public single-cell endometriosis studies, GSE213216
and GSE179640, as candidate extension cohorts for the donor-level stromal contrasts of this
project. Reads the public GEO SOFT family files, GEO supplementary file listings and the
HECA authors' public donor and library metadata (which curate these two studies' clinical
fields), joins them at donor level, applies the fixed eligibility rules below and writes the
sample manifests, the eutopic donor ledger, the lesion-pair table and a receipt. No expression
value is read and nothing is downloaded by this script. Standard library only; run with -I."""
import argparse,csv,datetime,gzip,hashlib,json,pathlib,re
from collections import Counter,defaultdict

RULES={"eutopic_scrna":"single-cell gene-expression library of eutopic endometrium (bulk, organoid and hashing libraries excluded)",
       "natural_cycle":"HECA-curated Hormonal treatment empty or '-', and Binary Stage in {Proliferative, Secretory} (Hormones, Menopause and Menstrual excluded; same rule family as v1)",
       "disease":"study condition (GSE179640) or HECA-curated endometrial pathology containing the E token / Endometriosis_stage > 0 (GSE213216)",
       "age_required":True,
       "minimum_per_disease_group_floor":5,"minimum_per_disease_group_preferred":8,
       "paired_lesion_design_minimum_pairs":5,
       "verdict_logic":"usable for the eutopic disease contrast only if, within one phase, both disease groups reach the floor among natural-cycle eutopic scRNA donors with age; lesion-context usable for a paired design only if matched eutopic+lesion scRNA pairs reach the pair minimum"}
URLS={"esummary":"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=200213216,200179640&retmode=json",
      "GSE213216":{"soft":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213216/soft/GSE213216_family.soft.gz","filelist":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213216/suppl/filelist.txt"},
      "GSE179640":{"soft":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179640/soft/GSE179640_family.soft.gz","filelist":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179640/suppl/filelist.txt"}}

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()
def write_tsv(rows,path,cols):
    with open(path,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter="\t",extrasaction="ignore",lineterminator="\n");w.writeheader()
        for r in rows:w.writerow(r)
def read_csv(p):
    with open(p,newline="",encoding="utf-8",errors="replace") as f:return list(csv.DictReader(f))

def parse_soft(path):
    series={};samples={};cur=None;kind=None
    with gzip.open(path,"rt",encoding="utf-8",errors="replace") as f:
        for line in f:
            line=line.rstrip("\n")
            if line.startswith("^SERIES"):kind="series";continue
            if line.startswith("^SAMPLE"):kind="sample";cur=samples.setdefault(line.split("=")[1].strip(),defaultdict(list));continue
            if line.startswith("^PLATFORM"):kind="platform";cur=defaultdict(list);continue
            if line.startswith("!") and "=" in line:
                k,v=line[1:].split("=",1);k=k.strip();v=v.strip()
                if kind=="series":series.setdefault(k,[]).append(v)
                elif cur is not None:cur[k].append(v)
    return series,samples
def chars(s):
    d={}
    for c in s.get("Sample_characteristics_ch1",[]):
        if ":" in c:k,v=c.split(":",1);d[k.strip().lower()]=v.strip()
    return d
def filelist(path):
    out=[]
    # columns: Archive, Name, Time, Size, Type; the row whose first field is the literal
    # "Archive" describes the RAW archive itself, the header row has Name in the second field
    for line in pathlib.Path(path).read_text().splitlines():
        p=line.split("\t")
        if len(p)<5 or p[1]=="Name":continue
        try:size=int(p[3])
        except ValueError:continue
        out.append({"name":p[1],"bytes":size,"type":"archive" if p[0]=="Archive" else p[4]})
    return out

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root",type=pathlib.Path,default=pathlib.Path(__file__).resolve().parents[1])
    ap.add_argument("--metadata",default="data/metadata/cohort_extension",help="directory holding the fetched GEO metadata files")
    ap.add_argument("--heca",default="data/metadata/heca_author")
    ap.add_argument("--retrieved-utc",required=True,help="when the GEO metadata files were fetched")
    ap.add_argument("--git-head",required=True);a=ap.parse_args()
    root=a.root.resolve();M=root/a.metadata;H=root/a.heca
    inputs={}
    for p in [M/"esummary_gds.json",M/"GSE213216/GSE213216_family.soft.gz",M/"GSE213216/filelist.txt",M/"GSE179640/GSE179640_family.soft.gz",M/"GSE179640/filelist.txt",
              H/"metadata-donors.csv",H/"metadata-libraries.csv",H/"metadata-libraries_deconvoluted.csv"]:
        inputs[str(p.relative_to(root)).replace("\\","/")]={"sha256":sha(p),"bytes":p.stat().st_size}
    summ=json.loads((M/"esummary_gds.json").read_text())["result"]
    donors=read_csv(H/"metadata-donors.csv");libs=read_csv(H/"metadata-libraries.csv");dec=read_csv(H/"metadata-libraries_deconvoluted.csv")
    donor_by_id={d["Donor"]:d for d in donors}
    dec_by_gsm={};dec_by_sample={}
    for r in dec:
        if r.get("StudyName")=="Fonseca":dec_by_gsm[r["sample"]]=r
        if r.get("StudyName")=="Tan":dec_by_sample[r["sample"]]=r
    manifests={};series_info={};ledger=[];pairs={};verdicts={}
    for gse in ("GSE213216","GSE179640"):
        series,samples=parse_soft(M/gse/f"{gse}_family.soft.gz");fl=filelist(M/gse/"filelist.txt")
        rows=[]
        for gsm,s in samples.items():
            c=chars(s)
            r={"gsm":gsm,"title":s.get("Sample_title",[""])[0],"source_name":s.get("Sample_source_name_ch1",[""])[0],
               "library_strategy":s.get("Sample_library_strategy",[""])[0],"instrument":s.get("Sample_instrument_model",[""])[0],
               "supplementary_files":";".join(x.split("/")[-1] for k,v in s.items() if k.startswith("Sample_supplementary_file") for x in v)}
            for k,v in c.items():r["char_"+re.sub(r"\W+","_",k)]=v
            rows.append(r)
        rows.sort(key=lambda r:r["gsm"])
        cols=["gsm","title","source_name","library_strategy","instrument"]+sorted({k for r in rows for k in r if k.startswith("char_")})+["supplementary_files"]
        manifests[gse]={"rows":rows,"cols":cols}
        series_info[gse]={"title":series.get("Series_title",[""])[0],"pubmed_ids":series.get("Series_pubmed_id",[]),"platforms":series.get("Series_platform_id",[]),
                          "submission_date":series.get("Series_submission_date",[""])[0],"last_update_date":series.get("Series_last_update_date",[""])[0],
                          "n_samples":len(samples),"summary":" ".join(series.get("Series_summary",[]))[:1500],"overall_design":" ".join(series.get("Series_overall_design",[]))[:1500],
                          "esummary":{k:summ.get("200"+gse[3:],{}).get(k) for k in ("n_samples","pubmedids","gdstype","summary","ftplink","suppfile")},
                          "supplementary_bytes_total":sum(f["bytes"] for f in fl if f["type"]!="archive"),"supplementary_files":fl,
                          "characteristic_keys":dict(Counter(k for r in rows for k in r if k.startswith("char_")))}
        # ---- eutopic single-cell donors with the HECA-curated clinical fields ----
        if gse=="GSE213216":
            eut=[r for r in rows if r.get("char_tissue")=="Eutopic Endometrium" and r["library_strategy"]=="RNA-Seq"]
            for r in eut:
                d=dec_by_gsm.get(r["gsm"],{});don=donor_by_id.get(d.get("genotype",""),{})
                path=don.get("Endometrial pathology_A-F-H-P","");stage=don.get("Endometriosis_stage","")
                disease="ENDO" if ("E" in [t.strip() for t in path.split("+")] or (stage.isdigit() and int(stage)>0)) else ("non-ENDO" if don else "")
                ledger.append({"study":gse,"gsm":r["gsm"],"geo_title":r["title"],"geo_patient_id":r.get("char_patient_id",""),"heca_library":d.get("Library_genotype",""),
                               "heca_donor":d.get("genotype",""),"disease":disease,"disease_source":"HECA-curated endometrial pathology / Endometriosis_stage","binary_stage":don.get("Binary Stage",""),
                               "stage":don.get("Stage",""),"hormonal_treatment":don.get("Hormonal treatment",""),"age":don.get("Age",""),"endometrial_pathology":path,
                               "endometriosis_stage":stage,"heca_biopsy_type":d.get("Biopsy_type",""),"heca_discard_lowqc":d.get("discard_lowQC",""),"heca_qc":d.get("QC",""),
                               "supplementary_bytes":next((f["bytes"] for f in fl if f["name"].startswith(r["gsm"]+"_")),None)})
            lesion_gsm=[r for r in rows if r.get("char_tissue") in ("Endometriosis","Endometrioma") and r["library_strategy"]=="RNA-Seq"]
            by_patient=defaultdict(lambda:{"eutopic":[],"lesion":[]})
            for r in eut:by_patient[r.get("char_patient_id","")]["eutopic"].append(r["gsm"])
            for r in lesion_gsm:by_patient[r.get("char_patient_id","")]["lesion"].append(r["gsm"]+":"+r.get("char_tissue","")+"/"+r.get("char_anatomic_location",""))
        else:
            eut=[r for r in rows if r.get("char_tissue")=="Endometrium" and r.get("char_method")=="scRNA-seq" and r.get("char_library_type")=="Gene Expression"]
            for r in eut:
                did=r["source_name"];don=donor_by_id.get(did,{});sample=r["supplementary_files"].split("_filtered")[0].split("_",1)[-1] if r["supplementary_files"] else ""
                d=dec_by_sample.get(sample,{})
                ledger.append({"study":gse,"gsm":r["gsm"],"geo_title":r["title"],"geo_patient_id":did,"heca_library":d.get("Library_genotype",""),"heca_donor":did if don else "",
                               "disease":"ENDO" if r.get("char_condition")=="Endometriosis" else ("non-ENDO" if r.get("char_condition")=="Control" else ""),"disease_source":"GEO condition",
                               "binary_stage":don.get("Binary Stage",""),"stage":don.get("Stage",""),"hormonal_treatment":don.get("Hormonal treatment",""),"age":don.get("Age",""),
                               "endometrial_pathology":don.get("Endometrial pathology_A-F-H-P",""),"endometriosis_stage":don.get("Endometriosis_stage",""),
                               "heca_biopsy_type":d.get("Biopsy_type",""),"heca_discard_lowqc":d.get("discard_lowQC",""),"heca_qc":d.get("QC",""),
                               "supplementary_bytes":next((f["bytes"] for f in fl if f["name"].startswith(r["gsm"]+"_")),None)})
            lesion_gsm=[r for r in rows if r.get("char_method")=="scRNA-seq" and r.get("char_library_type")=="Gene Expression" and r.get("char_sample_location") in ("Ectopic","Ectopic Ovary")]
            by_patient=defaultdict(lambda:{"eutopic":[],"lesion":[]})
            for r in eut:by_patient[r["source_name"]]["eutopic"].append(r["gsm"])
            for r in lesion_gsm:by_patient[r["source_name"]]["lesion"].append(r["gsm"]+":"+r.get("char_sample_location","")+"/"+r.get("char_tissue",""))
        pairs[gse]=[{"study":gse,"patient":p,"eutopic_scrna":";".join(v["eutopic"]),"lesion_scrna":";".join(v["lesion"]),"matched_pair":bool(v["eutopic"] and v["lesion"])} for p,v in sorted(by_patient.items())]
    # ---- eligibility and verdicts ----
    for r in ledger:
        horm=(r["hormonal_treatment"] or "").strip();natural=horm in ("","-") and r["binary_stage"] in ("Proliferative","Secretory")
        r["natural_cycle"]=natural;r["phase"]=r["binary_stage"] if natural else ""
        reasons=[]
        if not r["heca_donor"]:reasons.append("no curated clinical record")
        if horm not in ("","-"):reasons.append("hormonal treatment")
        if r["binary_stage"] not in ("Proliferative","Secretory"):reasons.append("stage "+(r["binary_stage"] or "missing"))
        if not r["age"]:reasons.append("age missing")
        if not r["disease"]:reasons.append("disease unresolved")
        r["eligible_eutopic_contrast"]=not reasons;r["exclusion_reasons"]="; ".join(reasons)
    for gse in ("GSE213216","GSE179640"):
        L=[r for r in ledger if r["study"]==gse]
        strata=Counter((r["phase"],r["disease"]) for r in L if r["eligible_eutopic_contrast"])
        counts={f"{ph} {dz}":n for (ph,dz),n in sorted(strata.items())}
        usable=any(strata.get((ph,"ENDO"),0)>=RULES["minimum_per_disease_group_floor"] and strata.get((ph,"non-ENDO"),0)>=RULES["minimum_per_disease_group_floor"] for ph in ("Proliferative","Secretory"))
        npairs=sum(1 for p in pairs[gse] if p["matched_pair"])
        nat_pairs=sum(1 for p in pairs[gse] if p["matched_pair"] and any(r["natural_cycle"] for r in L if r["geo_patient_id"]==p["patient"]))
        verdicts[gse]={"eutopic_scrna_donors":len(L),"eligible_natural_cycle_eutopic_donors_by_phase_and_disease":counts,
                       "exclusions":dict(Counter(r["exclusion_reasons"] for r in L if not r["eligible_eutopic_contrast"])),
                       "usable_for_eutopic_disease_contrast":usable,"matched_eutopic_lesion_scrna_pairs":npairs,"matched_pairs_with_natural_cycle_eutopic":nat_pairs,
                       "usable_for_paired_lesion_design_any_treatment":npairs>=RULES["paired_lesion_design_minimum_pairs"],
                       "usable_for_paired_lesion_design_natural_cycle":nat_pairs>=RULES["paired_lesion_design_minimum_pairs"]}
    heca_overlap={"heca_library_rows_by_study":dict(Counter(r.get("StudyName") for r in libs)),
                  "note":"both studies are external datasets integrated into the HECA cells atlas; the v1 and positive-control HECA cohorts use only the newly generated nuclei (E-MTAB-14039), so no donor of either study is a v1 HECA donor; HECA label vocabularies were derived on an atlas that includes them",
                  "fonseca_libraries_analysed_flag":dict(Counter(r.get("analysed") for r in libs if r.get("StudyName")=="Fonseca")),
                  "tan_libraries_analysed_flag":dict(Counter(r.get("analysed") for r in libs if r.get("StudyName")=="Tan"))}
    out_dir=root/"manifests";out_dir.mkdir(exist_ok=True);outputs={}
    for gse,m in manifests.items():
        p=out_dir/f"cohort_extension_{gse}_samples.tsv";write_tsv(m["rows"],p,m["cols"]);outputs[str(p.relative_to(root)).replace("\\","/")]=sha(p)
    lcols=["study","gsm","geo_title","geo_patient_id","heca_library","heca_donor","disease","disease_source","binary_stage","stage","phase","hormonal_treatment","age","endometrial_pathology","endometriosis_stage","natural_cycle","eligible_eutopic_contrast","exclusion_reasons","heca_biopsy_type","heca_discard_lowqc","heca_qc","supplementary_bytes"]
    p=out_dir/"cohort_extension_eutopic_donor_ledger.tsv";write_tsv(ledger,p,lcols);outputs[str(p.relative_to(root)).replace("\\","/")]=sha(p)
    p=out_dir/"cohort_extension_lesion_pairs.tsv";write_tsv([x for g in pairs.values() for x in g],p,["study","patient","eutopic_scrna","lesion_scrna","matched_pair"]);outputs[str(p.relative_to(root)).replace("\\","/")]=sha(p)
    receipt={"schema_version":1,"phase":"COHORT_EXTENSION_FEASIBILITY_METADATA_ONLY","created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
             "git_head":a.git_head,"retrieved_utc":a.retrieved_utc,"urls":URLS,"inputs":inputs,"rules":RULES,"series":series_info,"verdicts":verdicts,
             "heca_overlap":heca_overlap,"outputs_sha256":outputs,
             "attestations":{"expression_data_read":False,"expression_data_downloaded":False,"clinical_fields_source":"GEO sample characteristics and the HECA authors' public curation of these studies; not re-derived from the primary publications' supplements in this step"}}
    rp=root/"receipts/cohort_extension_metadata_audit.json";rp.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"verdicts":verdicts,"heca_overlap":heca_overlap["heca_library_rows_by_study"],"outputs":list(outputs)},indent=2))
if __name__=="__main__":main()
