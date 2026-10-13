#!/usr/bin/env python3
"""Bounded official source retrieval; preserve existing validated source files."""
import argparse,datetime,hashlib,json,os,shutil,urllib.request,zipfile
from pathlib import Path
SOURCES={
 'article':('https://pmc.ncbi.nlm.nih.gov/articles/PMC10042257/','data/genetics/pmc.txt',2*1024**2),
 'supplement':('https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-023-01323-z/MediaObjects/41588_2023_1323_MOESM4_ESM.xlsx','data/genetics/41588_2023_1323_MOESM4_ESM.xlsx',130*1024**2)}
def digest(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(1024**2):h.update(b)
 return h.hexdigest()
def validate(name,p):
 if name=='supplement':
  assert zipfile.is_zipfile(p),'Not XLSX'
  with zipfile.ZipFile(p) as z:assert 'xl/workbook.xml' in z.namelist()
 else:
  b=p.read_bytes();assert b'The genetic basis of endometriosis' in b and b'GCST90205183' in b,'Not primary article'
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--fetch-missing',action='store_true');args=parser.parse_args()
 manifest=Path('manifests/genetic_sources.tsv');expected={}
 if manifest.exists():
  import csv
  expected={r['path']:r['sha256'] for r in csv.DictReader(manifest.open(),delimiter='\t') if r['status']=='accepted'}
 rows=[]
 for name,(url,path,limit) in SOURCES.items():
  p=Path(path);reused=p.exists()
  if not reused:
   assert args.fetch_missing,'Missing source; use --fetch-missing'
   assert shutil.disk_usage('.').free-limit>20*1024**3,'Disk reserve would fail'
   p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.part');assert not tmp.exists(),'Partial already exists; inspect it'
   n=0
   with urllib.request.urlopen(url,timeout=90) as r,tmp.open('wb') as f:
    while b:=r.read(1024**2):
     n+=len(b);assert n<=limit,'Download cap exceeded';f.write(b)
   validate(name,tmp);os.replace(tmp,p)
  validate(name,p);h=digest(p)
  if path in expected:assert h==expected[path],'Source hash changed; preserve and investigate'
  rows.append(dict(name=name,path=path,url=url,bytes=p.stat().st_size,sha256=h,reused=reused))
 print(json.dumps(dict(checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sources=rows),indent=2))
if __name__=='__main__':main()
