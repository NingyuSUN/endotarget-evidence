"""Fetch only the fixed official HGNC object; existing bytes must match the pin."""
import argparse,base64,csv,hashlib,io,json,pathlib,urllib.request
ROOT=pathlib.Path(__file__).resolve().parents[1]
def validate_annotation(blob,pin):
    if len(blob)!=pin['bytes'] or hashlib.sha256(blob).hexdigest()!=pin['sha256']:
        raise ValueError('HGNC size or SHA-256 differs from frozen annotation')
    if base64.b64encode(hashlib.md5(blob).digest()).decode()!=pin['official_GCS_md5']:
        raise ValueError('HGNC official object MD5 mismatch')
    rows=list(csv.DictReader(io.StringIO(blob.decode('utf-8-sig')),delimiter='\t'))
    required={'hgnc_id','symbol','status','ensembl_gene_id','prev_symbol'}
    if not rows or not required<=set(rows[0]) or any(r['status']!='Approved' for r in rows):
        raise ValueError('HGNC approved annotation schema failed')
    return len(rows)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--reference',default='config/gene_mapping_reference.json');a=ap.parse_args()
    pin=json.loads((ROOT/a.reference).read_text());dest=ROOT/pin['path'];dest.parent.mkdir(parents=True,exist_ok=True)
    reused=dest.exists()
    if reused:blob=dest.read_bytes()
    else:
        with urllib.request.urlopen(pin['source_url'],timeout=60) as response:
            blob=response.read(25*1024*1024+1)
        if len(blob)>25*1024*1024:raise ValueError('HGNC exceeds 25MiB retrieval limit')
    n=validate_annotation(blob,pin)
    if not reused:
        tmp=dest.with_suffix(dest.suffix+'.partial');tmp.write_bytes(blob);tmp.replace(dest)
    report={**pin,'verified_rows':n,'reused_matching_pinned_bytes':reused,'validation':'fixed generation URL, bytes, official MD5, frozen SHA256 and approved schema'}
    (ROOT/'receipts/gene_annotation_hgnc_replay.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
