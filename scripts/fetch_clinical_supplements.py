#!/usr/bin/env python3
"""Bounded official supplement retrieval; no disease-result sheet inspection."""
import hashlib,json,signal,time,urllib.request,zipfile,io
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
CAP=128*1024*1024
attempts=[]
class OverallDeadline(RuntimeError): pass
def deadline(*args): raise OverallDeadline("overall deadline 180s exceeded")
signal.signal(signal.SIGALRM,deadline)
def fetch(url,target,package=False):
 a={"url":url,"target":target,"attempt_utc":datetime.now(timezone.utc).isoformat(),"cap_bytes":CAP,"deadline_seconds":180,"bytes_received":0}
 start=time.monotonic(); signal.alarm(180)
 try:
  req=urllib.request.Request(url,headers={"User-Agent":"EndoTarget public metadata audit"})
  with urllib.request.urlopen(req,timeout=40) as r:
   a.update(final_url=r.url,content_type=r.headers.get("Content-Type"),content_length=r.headers.get("Content-Length"))
   blocks=[]
   while True:
    b=r.read(256*1024)
    if not b:break
    a["bytes_received"]+=len(b)
    if a["bytes_received"]>CAP:raise ValueError("download exceeds 128MiB cap")
    blocks.append(b)
  raw=b"".join(blocks)
  if not raw.startswith(b"PK"):raise ValueError("response is not ZIP/XLSX; HTML/CAPTCHA is not bypassed")
  if package:
   with zipfile.ZipFile(io.BytesIO(raw)) as z:
    matches=[n for n in z.namelist() if Path(n).name==target]
    a["package_members"]=len(z.namelist())
    if len(matches)!=1:raise ValueError("exact target not unique or absent in official package")
    member=z.getinfo(matches[0])
    if member.file_size>CAP:raise ValueError("selected workbook exceeds cap")
    raw=z.read(matches[0]);a["selected_member"]=matches[0]
  with zipfile.ZipFile(io.BytesIO(raw)) as z:
   for n in ["[Content_Types].xml","xl/workbook.xml","xl/_rels/workbook.xml.rels"]:
    if n not in z.namelist():raise ValueError("missing XLSX component: "+n)
   bad=z.testzip()
   if bad:raise ValueError("CRC failure: "+bad)
  out=ROOT/"data/clinical"/target
  out.write_bytes(raw)
  a.update(status="retrieved",saved_path=str(out.relative_to(ROOT)),workbook_bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
  return True
 except Exception as e:a.update(status="failed",error=type(e).__name__+": "+str(e));return False
 finally:
  signal.alarm(0);a["elapsed_seconds"]=round(time.monotonic()-start,3);attempts.append(a)
  (ROOT/"receipts/clinical_recovery_downloads.json").write_text(json.dumps(attempts,indent=2)+"\n")
  print(json.dumps(a),flush=True)
def main():
 for pmc,target in [("PMC13385355","41467_2026_73020_MOESM3_ESM.xlsx"),("PMC11387200","41588_2024_1873_MOESM3_ESM.xlsx")]:
  if fetch(f"https://pmc.ncbi.nlm.nih.gov/articles/{pmc}/bin/{target}",target):continue
  fetch(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/supplementaryFiles",target,True)

def fetch_streamed_package(pmc,target):
 """Read local ZIP records; stop immediately after the exact selected workbook."""
 import struct,zlib
 url=f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/supplementaryFiles"
 a={"url":url,"target":target,"method":"sequential ZIP local records, stop at target","attempt_utc":datetime.now(timezone.utc).isoformat(),"cap_bytes":CAP,"deadline_seconds":180,"bytes_received":0,"visited_members":[]}
 start=time.monotonic();signal.alarm(180)
 try:
  with urllib.request.urlopen(urllib.request.Request(url,headers={"User-Agent":"EndoTarget public metadata audit"}),timeout=40) as r:
   a.update(final_url=r.url,content_type=r.headers.get("Content-Type"))
   buf=b""
   def take(n):
    nonlocal buf
    while len(buf)<n:
     b=r.read(16384);a["bytes_received"]+=len(b)
     if not b:raise ValueError("premature ZIP end")
     if a["bytes_received"]>CAP:raise ValueError("package cap exceeded")
     buf+=b
    b,buf=buf[:n],buf[n:];return b
   while True:
    if take(4)!=b"PK\x03\x04":raise ValueError("target absent or unsupported ZIP record")
    h=struct.unpack("<5H3I2H",take(26));flags,method,csize,usize,nlen,xlen=h[1],h[2],h[6],h[7],h[8],h[9]
    name=take(nlen).decode("utf-8");take(xlen);a["visited_members"].append(name)
    selected=Path(name).name==target;raw=[];total=0
    if method==8:
     d=zlib.decompressobj(-15)
     remaining=csize if not flags&8 else None
     while not d.eof:
      count=min(16384,remaining) if remaining is not None else 16384
      if count<=0:raise ValueError("invalid compressed ZIP length")
      block=take(count);decoded=d.decompress(block);total+=len(decoded)
      if selected:
       if total>CAP:raise ValueError("workbook cap exceeded")
       raw.append(decoded)
      if remaining is not None:remaining-=len(block)
      if d.eof:buf=d.unused_data+buf
    elif method==0 and not flags&8:
     decoded=take(csize);total=len(decoded)
     if selected:raw.append(decoded)
    else:raise ValueError("unsupported streaming ZIP compression/descriptor")
    if flags&8:
     first=take(4)
     take(12 if first==b"PK\x07\x08" else 8)
    if selected:
     data=b"".join(raw)
     with zipfile.ZipFile(io.BytesIO(data)) as z:
      if z.testzip():raise ValueError("workbook CRC failure")
      if not {"xl/workbook.xml","xl/_rels/workbook.xml.rels","[Content_Types].xml"}.issubset(z.namelist()):raise ValueError("XLSX components absent")
     out=ROOT/"data/clinical"/target;out.write_bytes(data)
     a.update(status="retrieved",sha256=hashlib.sha256(data).hexdigest(),workbook_bytes=len(data),saved_path=str(out.relative_to(ROOT)))
     break
 except Exception as e:a.update(status="failed",error=type(e).__name__+": "+str(e))
 finally:
  signal.alarm(0);a["elapsed_seconds"]=round(time.monotonic()-start,3)
  p=ROOT/"receipts/clinical_recovery_stream.json"
  old=json.loads(p.read_text()) if p.exists() else [];old.append(a);p.write_text(json.dumps(old,indent=2)+"\n")
  print(json.dumps(a),flush=True)

if __name__=="__main__":
 import argparse
 parser=argparse.ArgumentParser()
 parser.add_argument("--stream",choices=["discovery","heca"])
 args=parser.parse_args()
 targets={"discovery":("PMC13385355","41467_2026_73020_MOESM3_ESM.xlsx"),"heca":("PMC11387200","41588_2024_1873_MOESM3_ESM.xlsx")}
 if args.stream:fetch_streamed_package(*targets[args.stream])
 else:main()
