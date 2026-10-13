#!/usr/bin/env python3
"""Bounded public metadata retrieval and conservative donor audit; Python >=3.10."""
import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
import os
import re
import ssl
import signal
from contextlib import contextmanager, nullcontext
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

MAX_BYTES = 10 * 1024 * 1024
TIMEOUT = 30
TOTAL_TIMEOUT = 120
DONOR = "Characteristics[individual]"
ASSAY = "Characteristics[10x technology]"
DISEASE = "Characteristics[disease]"
PHASE = "Characteristics[menstrual_cycle_stage]"
SOURCES = [
    {"id":"sdrf", "url":"https://www.ebi.ac.uk/biostudies/files/E-MTAB-14039/E-MTAB-14039.sdrf.txt", "kind":"sdrf", "filename":"E-MTAB-14039.sdrf.txt", "provenance":"manifests/source_inventory.csv"},
    {"id":"biostudies", "url":"https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-14039", "kind":"biostudies_json", "filename":"E-MTAB-14039.json", "provenance":"docs/DATA_SOURCES.md"},
    {"id":"geo_soft", "url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE266nnn/GSE266265/soft/GSE266265_family.soft.gz", "kind":"geo_soft_gz", "filename":"GSE266265_family.soft.gz", "provenance":"GEO official accession FTP family directory; alternate after documented GEO record CAPTCHA rejection"},
]

def utc():
    return datetime.now(timezone.utc).isoformat()

def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    with part.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(part, path)

def parse_sdrf(data):
    text = data.decode("utf-8-sig")
    if text.lstrip().startswith("<") or "\x00" in text:
        raise ValueError("SDRF payload is HTML/XML or contains NUL")
    rows = list(csv.reader(io.StringIO(text), delimiter="\t", strict=True))
    if len(rows) < 2 or not rows[0] or rows[0][0] != "Source Name":
        raise ValueError("SDRF needs Source Name header and data rows")
    fields = rows[0]
    required = {DONOR, ASSAY, DISEASE, PHASE}
    if not required.issubset(fields):
        raise ValueError("Missing required SDRF columns: " + repr(sorted(required-set(fields))) + "; observed=" + repr(fields))
    if any(fields.count(field) != 1 for field in required):
        raise ValueError("Duplicate required SDRF column")
    if any(len(row) != len(fields) for row in rows[1:]):
        raise ValueError("Ragged SDRF rows")
    # SDRF legitimately repeats ontology/protocol headings; retain ordered columns.
    return fields, rows[1:]

def validate_payload(data, kind):
    if not data or len(data) > MAX_BYTES:
        raise ValueError("Empty or over-limit metadata")
    if kind == "sdrf":
        fields, rows = parse_sdrf(data)
        return {"columns":len(fields), "rows":len(rows)}
    if kind == "biostudies_json":
        try:
            obj = json.loads(data)
        except (ValueError, UnicodeError) as exc:
            raise ValueError("Invalid BioStudies JSON") from exc
        if not isinstance(obj, dict) or obj.get("accno") != "E-MTAB-14039" or not isinstance(obj.get("section"),dict):
            raise ValueError("Wrong BioStudies accession/schema")
        return {"accession":obj["accno"], "section_type":obj["section"].get("type")}
    if kind == "geo_soft_gz":
        try:
            with gzip.GzipFile(fileobj=io.BytesIO(data)) as handle:
                raw = handle.read(MAX_BYTES+1)
            if len(raw)>MAX_BYTES:
                raise ValueError("Decompressed SOFT exceeds byte limit")
            text = raw.decode("utf-8")
        except (OSError, EOFError, UnicodeError) as exc:
            raise ValueError("Invalid gzip SOFT") from exc
        if not re.search(r"^\^SERIES = GSE266265$",text,re.M) or "!Series_title = " not in text:
            raise ValueError("Wrong GEO SOFT accession/schema")
        samples = re.findall(r"^\^SAMPLE = (GSM[0-9]+)$",text,re.M)
        if not samples:
            raise ValueError("Family SOFT contains no sample records")
        return {"accession":"GSE266265", "decompressed_bytes":len(raw),
                "sample_records":len(samples),"sample_accessions":samples,
                "note":"GEO sample records are not verified eligible donors"}
    if kind == "geo_record":
        text = data.decode("utf-8")
        if not ("GSE266265" in text and re.search(r"Series\s+GSE266265",text) and "GSM" in text and "Status" in text):
            raise ValueError("Not the expected GEO study record")
        return {"accession":"GSE266265", "representation":"official HTML record; no per-donor eligibility inferred"}
    raise ValueError("Unknown metadata type")

@contextmanager
def overall_deadline(seconds):
    """Interrupt blocking open/read on the SRV01 Unix main-thread CLI."""
    if signal.getitimer(signal.ITIMER_REAL)[0]:
        raise RuntimeError("Refusing to replace an existing process alarm")
    previous_handler = signal.getsignal(signal.SIGALRM)
    def expired(signum, frame):
        raise TimeoutError("Download exceeded overall deadline")
    signal.signal(signal.SIGALRM, expired)
    try:
        signal.setitimer(signal.ITIMER_REAL, seconds)
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)

def stream_to_part(response, part, limit=MAX_BYTES, *, enforce_deadline=True):
    if not 0 < limit <= MAX_BYTES:
        raise ValueError("Invalid streaming limit")
    size = 0
    digest = hashlib.sha256()
    # fetch_one supplies its own deadline covering both open and read.
    guard = overall_deadline(TOTAL_TIMEOUT) if enforce_deadline else nullcontext()
    # x mode leaves preexisting partial files untouched.
    with part.open("xb") as handle:
        try:
            with guard:
                while True:
                    chunk = response.read(min(65536, limit-size+1))
                    if not chunk:
                        break
                    size += len(chunk)
                    if size > limit:
                        raise ValueError("Download exceeds byte limit")
                    handle.write(chunk)
                    digest.update(chunk)
                handle.flush()
                os.fsync(handle.fileno())
        except BaseException:
            handle.close()
            part.unlink()
            raise
    return size, digest.hexdigest()

class HTTPSOnlyRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if urllib.parse.urlsplit(newurl).scheme != "https":
            raise ValueError("Refusing non-HTTPS redirect")
        return super().redirect_request(req, fp, code, msg, headers, newurl)

def fetch_one(source, root, refresh=False):
    url = source["url"]
    if urllib.parse.urlsplit(url).scheme != "https":
        raise ValueError("Only HTTPS sources allowed")
    target = root/"data/metadata"/source["filename"]
    receipt_path = root/"receipts"/("m0_metadata_"+source["id"]+".json")
    target.parent.mkdir(parents=True, exist_ok=True)
    old = None
    if target.exists():
        if not receipt_path.exists():
            raise ValueError("Existing file has no receipt; preserved, not overwritten")
        old = json.loads(receipt_path.read_text())
        if target.stat().st_size > MAX_BYTES:
            raise ValueError("Existing metadata exceeds byte limit; preserved")
        data = target.read_bytes()
        if (hashlib.sha256(data).hexdigest() != old.get("sha256") or len(data) != old.get("bytes")
                or old.get("url") != url or old.get("kind") != source["kind"]):
            raise ValueError("Existing metadata does not match receipt; preserved")
        validate_payload(data, source["kind"])
        if not refresh:
            return dict(old, status="reused_validated", cache_checked_utc=utc())
    part = target.with_name(target.name+".part")
    opener = urllib.request.build_opener(
        urllib.request.HTTPSHandler(context=ssl.create_default_context()), HTTPSOnlyRedirect())
    request = urllib.request.Request(url, headers={"User-Agent":"EndoTarget-MetadataAudit/0.1", "Accept-Encoding":"identity"})
    try:
        with overall_deadline(TOTAL_TIMEOUT):
            with opener.open(request, timeout=TIMEOUT) as response:
                final_url = response.geturl()
                if urllib.parse.urlsplit(final_url).scheme != "https":
                    raise ValueError("Final URL is not HTTPS")
                length = response.headers.get("Content-Length")
                if length is not None and int(length) > MAX_BYTES:
                    raise ValueError("Content-Length exceeds byte limit")
                size, sha = stream_to_part(response, part, enforce_deadline=False)
                content_type = response.headers.get("Content-Type")
        schema = validate_payload(part.read_bytes(), source["kind"])
    except BaseException:
        # Do not delete a stale .part we did not create.
        if "size" in locals() and part.exists():
            rejected = target.with_name(target.name+".rejected."+sha)
            if rejected.exists():
                if hashlib.sha256(rejected.read_bytes()).hexdigest() != sha:
                    raise ValueError("Existing rejected-payload hash mismatch")
                part.unlink()
            else:
                os.replace(part, rejected)
            atomic_json(root/"receipts"/("m0_metadata_"+source["id"]+".rejected."+sha+".json"),
                        dict(source,status="SCHEMA_REJECTED",path=str(rejected),sha256=sha,bytes=size,
                             retrieved_utc=utc(),final_url=final_url))
        raise
    receipt = dict(source, url=url, final_url=final_url, bytes=size, sha256=sha,
                   retrieved_utc=utc(), content_type=content_type, schema=schema,
                   https_certificate_verification=True, byte_limit=MAX_BYTES,
                   socket_timeout_seconds=TIMEOUT, overall_deadline_seconds=TOTAL_TIMEOUT,
                   path=str(target), status="downloaded")
    if old is not None:
        if sha == old["sha256"]:
            part.unlink()
            return dict(old, status="remote_unchanged", remote_checked_utc=receipt["retrieved_utc"])
        changed = target.with_name(target.name+".changed."+sha)
        if changed.exists():
            if hashlib.sha256(changed.read_bytes()).hexdigest() != sha:
                raise ValueError("Changed-file path contains unexpected bytes")
            part.unlink()
        else:
            os.replace(part, changed)
        receipt.update(status="source_changed", changed_path=str(changed), original_sha256=old["sha256"])
        atomic_json(root/"receipts"/("m0_metadata_"+source["id"]+".changed."+sha+".json"),receipt)
        return receipt
    os.replace(part, target)
    atomic_json(receipt_path,receipt)
    return receipt

def audit_sdrf(data):
    fields, rows = parse_sdrf(data)
    indices = {field:fields.index(field) for field in [DONOR,ASSAY,DISEASE,PHASE]}
    clinical = [(i,f) for i,f in enumerate(fields) if f.startswith("Characteristics[") and f not in (DONOR, ASSAY, "Characteristics[library_id]")]
    groups = collections.defaultdict(list)
    missing_identity = []
    for number,row in enumerate(rows,2):
        key = (row[indices[ASSAY]],row[indices[DONOR]])
        if not all(v.strip() for v in key):
            missing_identity.append({"sdrf_line":number, "assay":key[0], "declared_donor":key[1]})
        else:
            groups[key].append((number,row))
    donors, conflicts, missing, naming = [], [], [], []
    strata, candidates, assay_counts = collections.Counter(),collections.Counter(),collections.Counter()
    for (assay,donor), members in sorted(groups.items()):
        values = {f:sorted({r[i] for _,r in members}) for i,f in clinical}
        bad = {f:v for f,v in values.items() if len(v)>1}
        empty = {f:v for f,v in values.items() if any(vv.strip().lower() in ("","unknown","not available","not provided","na","n/a") for vv in v)}
        entry = {"assay":assay,"declared_donor":donor,"library_rows":len(members),
                 "sdrf_lines":[n for n,_ in members], "source_values":values, "library_ids":sorted({r[fields.index("Characteristics[library_id]")] for _,r in members}) if "Characteristics[library_id]" in fields else []}
        donors.append(entry)
        assay_counts[assay] += 1
        if bad:
            conflicts.append({"assay":assay,"declared_donor":donor,"fields":bad,"sdrf_lines":entry["sdrf_lines"]})
        if empty:
            missing.append({"assay":assay,"declared_donor":donor,"fields":empty})
        diseases, phases = values[DISEASE], values[PHASE]
        if len(diseases)==len(phases)==1:
            strata[(assay,diseases[0],phases[0])] += 1
        # Candidate means only exact source phase+disease match, not eligibility.
        if not bad and phases==["Proliferative"] and diseases[0] in ("endometriosis","normal"):
            candidates[(assay,diseases[0])] += 1
        for line,row in members:
            for index,field in enumerate(fields):
                if "file" in field.lower() and row[index]:
                    # Evidence only: report declared ID absent from filename.
                    filename = urllib.parse.unquote(row[index].rsplit("/",1)[-1])
                    if not re.search(r"(?<![A-Za-z0-9])"+re.escape(donor)+r"(?![A-Za-z0-9])",filename):
                        naming.append({"sdrf_line":line,"assay":assay,"declared_donor":donor,
                                       "field":field,"source_value":row[index],
                                       "issue":"declared_ID_not_literal_filename_token",
                                       "status":"UNRESOLVED_NAMING_EVIDENCE","reassigned":False})
    readiness = {
        "stage0_passed":False,"status":"OPEN","final_eligibility":"NOT_VERIFIED",
        "eligible_discovery_donors":None,"eligible_replication_donors":None,
        "metadata_only":True,
        "required_unverified_gates":["cell_to_donor_mapping","stromal_cells_per_donor",
            "hormone_exposure_and_clinical_exclusions","raw_count_layer",
            "age_batch_model_identifiability","independence_ledger","resource_feasibility","source_terms"],
        "fixed_thresholds":{"discovery_min_per_group":8,"replication_min_per_group":5,"stromal_cells_min":50},
        "metadata_conflict_units":len(conflicts), "missing_identity_rows":len(missing_identity),
        "candidate_counts_are_final_eligibility":False
    }
    return {"schema_version":1, "status":"METADATA_ONLY","field_names":fields,
            "deduplication_key":[ASSAY,DONOR],"library_rows":len(rows),"donor_units":len(donors),
            "donor_units_note":"Assay-specific declared IDs; do not sum assays as independent people.",
            "donors_by_assay":[{"assay":a,"donors":n} for a,n in sorted(assay_counts.items())],
            "donor_strata":[{"assay":a,"disease":d,"phase":p,"donors":n} for (a,d,p),n in sorted(strata.items())],
            "candidate_proliferative_counts":[{"assay":a,"disease":d,"donors":n} for (a,d),n in sorted(candidates.items())],
            "candidate_definition":"Exact Proliferative and endometriosis/normal labels; no within-donor clinical conflicts. Hormone, cell/QC and count eligibility NOT VERIFIED.",
            "within_donor_conflicts":conflicts,"missing_clinical_fields":missing,
            "rows_missing_identity":missing_identity,"filename_naming_evidence":naming,
            "filename_check_limit":"Literal-token absence is a broad screen, not a donor inference; opaque filenames can trigger it.",
            "donors":donors, "readiness":readiness}

def run(root, refresh=False):
    records = []
    for source in SOURCES:
        try:
            records.append(fetch_one(source,root,refresh))
        except Exception as exc:
            records.append(dict(source,status="FAILED",error_type=type(exc).__name__,error=str(exc),failed_utc=utc()))
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    run_receipt = root/"receipts"/("m0_metadata_run_"+stamp+".json")
    atomic_json(run_receipt,{"run_utc":utc(),"host_execution":"SRV01 /home/analyst/te","records":records})
    usable = [r for r in records if r["status"] not in ("FAILED","source_changed")]
    # Update lock only for a complete unchanged/initial set; failed refresh preserves old pin.
    if len(usable)==len(SOURCES):
        lock = root/"manifests/inputs.lock.tsv"
        part = lock.with_name(lock.name+".part")
        with part.open("x",newline="",encoding="utf-8") as handle:
            writer = csv.writer(handle,delimiter="\t",lineterminator="\n")
            writer.writerow(["source_id","path","kind","bytes","sha256","retrieved_utc","url","final_url"])
            for r in usable:
                writer.writerow([r["id"],r["path"],r["kind"],r["bytes"],r["sha256"],r["retrieved_utc"],r["url"],r["final_url"]])
        os.replace(part,lock)
    sdrf = next((r for r in usable if r["id"]=="sdrf"), None)
    if sdrf:
        audit = audit_sdrf(Path(sdrf["path"]).read_bytes())
        audit.update(audited_utc=utc(),source_sha256=sdrf["sha256"],source_bytes=sdrf["bytes"],
                     source_url=sdrf["url"],acquisition_complete=len(usable)==len(SOURCES))
        atomic_json(root/"reports/m0_metadata_audit.json",audit)
        atomic_json(root/"reports/m0_metadata_readiness.json",audit["readiness"])
        print(json.dumps({k:audit[k] for k in ["library_rows","donor_units","donors_by_assay","candidate_proliferative_counts","readiness"]},indent=2))
    print(json.dumps({"run_receipt":str(run_receipt),"downloads":[{"id":r["id"],"status":r["status"],"error":r.get("error")} for r in records]},indent=2))
    return 0 if len(usable)==len(SOURCES) else 1

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument("--refresh",action="store_true",help="Re-fetch, detect remote changes, and preserve original pins.")
    args = parser.parse_args()
    return run(args.root.resolve(),args.refresh)

if __name__=="__main__":
    raise SystemExit(main())
