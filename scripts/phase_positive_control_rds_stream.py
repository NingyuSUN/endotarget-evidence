"""Read-only XDR RDS structural extractor; project port of the archived author parser.

Parsing logic is identical to the archived `rds_stream.py` (SHA-256
4154b22c7921d1edb0b9c3e01663f1a6fc3ed48f720f88bf09ed0ddfbd67b3de, OSF route inputs).
Three operational fixes, documented in docs/PHASE_POSITIVE_CONTROL_DISCOVERY_COUNTS_PLAN.md:
1. the resident-memory guard converts Linux ru_maxrss (KiB) to bytes before comparing;
2. the memory ceiling is a parameter (default 24 GiB) instead of a fixed 35 GiB;
3. the output-volume disk reserve is a parameter (default 75 GiB).
Numeric payloads are copied in <=8 MiB chunks to disk or discarded. The parser never
evaluates R code or constructs Seurat objects; unsupported types abort rather than guess.
"""
import argparse, bz2, gzip, hashlib, json, lzma, pickle, platform, resource, shutil, struct, sys, time
from dataclasses import dataclass
from pathlib import Path
import numpy as np

ARCHIVED_SHA256='4154b22c7921d1edb0b9c3e01663f1a6fc3ed48f720f88bf09ed0ddfbd67b3de'
RSS_UNIT=1024 if sys.platform.startswith('linux') else 1  # ru_maxrss is KiB on Linux, bytes on macOS

@dataclass
class Node:
    type: int
    value: object = None
    attrs: object = None

@dataclass
class DiskVector:
    path: str
    dtype: str
    length: int
    sha256: str

def value(x):
    if isinstance(x,Node): return x.value
    return x

def attrs(x):
    return dict(x.attrs or []) if isinstance(x,Node) else {}

def rss_bytes():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*RSS_UNIT

class Reader:
    def __init__(self,stream,out,keep_all=False,rss_ceiling_bytes=24*1024**3,disk_reserve_bytes=75*1024**3):
        self.f=stream;self.out=Path(out);self.out.mkdir(parents=True,exist_ok=True)
        self.keep_all=keep_all; self.refs=[None];self.pos=0;self.nodes=0;self.payloads=[];self.start=time.time();self.last=0
        self.rss_ceiling=int(rss_ceiling_bytes);self.disk_reserve=int(disk_reserve_bytes)
    def read(self,n):
        if n<0 or n>8*1024*1024: raise ValueError(f'Unsafe read {n}')
        if n>1024*1024:self.guard()
        b=self.f.read(n)
        if len(b)!=n: raise EOFError(f'Expected {n} bytes at {self.pos}, got {len(b)}')
        self.pos+=n;return b
    def i(self): return struct.unpack('>i',self.read(4))[0]
    def length(self):
        n=self.i()
        if n==-1: n=(self.i()<<32)+(self.i()&0xffffffff)
        if n<0: raise ValueError('negative vector length')
        return n
    def guard(self):
        rss=rss_bytes()
        if rss>self.rss_ceiling: raise MemoryError(f'parser RSS {rss} bytes exceeds ceiling {self.rss_ceiling} bytes')
        if time.time()-self.last>30:
            self.last=time.time()
            print(json.dumps({'elapsed_s':round(time.time()-self.start),'stream_GB':round(self.pos/1e9,3),'peak_RSS_GiB':round(rss/1024**3,3),'nodes':self.nodes}),flush=True)
    def payload(self,t,path):
        n=self.length();dtype={10:'>i4',13:'>i4',14:'>f8',15:'>c16',24:'u1'}[t]; size=np.dtype(dtype).itemsize*n
        keep=self.keep_all or path.startswith('root/@/assays/[0]/@/counts') or path.startswith('root/@/meta.data')
        # Tiny dimensions/flags are needed even outside retained payloads.
        if size<=8192:
            return np.frombuffer(self.read(size),dtype=dtype).copy()
        h=hashlib.sha256(); target=None
        if keep:
            if shutil.disk_usage(self.out).free<size+self.disk_reserve: raise OSError(f'{self.disk_reserve/1024**3:g} GiB disk reserve')
            target=self.out/f'vector_{len(self.payloads):06d}.bin'
            writer=target.open('wb')
        try:
            remain=size
            while remain:
                b=self.read(min(remain,8*1024**2));remain-=len(b)
                if keep: writer.write(b);h.update(b)
                self.guard()
        finally:
            if keep: writer.close()
        rec={'path':path,'type':t,'length':n,'bytes':size,'retained':keep,'file':str(target) if keep else None,'sha256':h.hexdigest() if keep else None}
        self.payloads.append(rec)
        return DiskVector(str(target),dtype,n,h.hexdigest()) if keep else {'discarded_length':n,'dtype':dtype}
    def stringvec(self,path):
        names=self.i()
        if names!=0: raise ValueError('Unexpected named namespace string vector')
        return [self.obj(path+'/namespace') for _ in range(self.length())]
    def pair(self,flags,path):
        result=[];index=0
        while flags&255 in (2,3,5,6,17):
            a=self.obj(path+'/pairattr') if flags&512 else None
            tag=self.obj(path+'/tag') if flags&1024 else None
            key=value(tag) if tag is not None else str(index)
            car=self.obj(path+'/'+str(key))
            result.append((str(key),car));index+=1
            flags=self.i()
        tail=self.obj(path+'/tail',flags)
        if tail is not None: result.append(('__tail__',tail))
        return result
    def bc_lang(self,t,path,reps):
        if t==243:return reps[self.i()]
        if t in (244,6,2,240,239):
            pos=None
            if t==244:pos=self.i();t=self.i()
            node=Node(t)
            if pos is not None:reps[pos]=node
            a=self.obj(path+'/attr') if t in (240,239) else None
            tag=self.obj(path+'/tag')
            car=self.bc_lang(self.i(),path+'/car',reps)
            cdr=self.bc_lang(self.i(),path+'/cdr',reps)
            node.value=(tag,car,cdr);node.attrs=a;return node
        return self.obj(path)
    def bc(self,path,reps=None):
        if reps is None:reps=[None]*self.length()
        code=self.obj(path+'/code'); constants=[]
        for j in range(self.length()):
            t=self.i()
            constants.append(self.bc(path+'/bc',reps) if t==21 else self.bc_lang(t,path+'/const',reps))
        return (code,constants)
    def obj(self,path,flags=None):
        self.nodes+=1
        if self.nodes%1000==0:self.guard()
        if flags is None:flags=self.i()
        t=flags&255
        if t in (0,254):return None
        if t in (241,242,250,251,252,253):return Node(t)
        if t==255:
            idx=(flags&0xffffffff)>>8
            if not idx:idx=self.i()
            return self.refs[idx]
        if t==1:
            obj=Node(t,value(self.obj(path+'/symbol')));self.refs.append(obj);return obj
        if t in (247,248,249):
            obj=Node(t,self.stringvec(path));self.refs.append(obj);return obj
        if t==238:
            info=self.obj(path+'/altinfo');state=self.obj(path+'/altstate');a=self.obj(path+'/altattr')
            return Node(t,(info,state),a)
        if t==4:
            locked=self.i();obj=Node(t);self.refs.append(obj)
            obj.value=(locked,self.obj(path+'/enclosure'),self.obj(path+'/frame'),self.obj(path+'/hashtable'))
            obj.attrs=self.obj(path+'/envattr');return obj
        if t in (2,3,5,6,17):return self.pair(flags,path)
        obj=Node(t)
        if t in (22,23):
            self.refs.append(obj)
            if t==22:obj.value=(self.obj(path+'/protected'),self.obj(path+'/tag'))
        elif t in (7,8,9):
            n=self.i();obj.value=None if n==-1 else self.read(n).decode('utf-8',errors='surrogateescape')
        elif t in (10,13,14,15,24):obj.value=self.payload(t,path)
        elif t in (16,19,20):
            n=self.length()
            if n>5_000_000: raise MemoryError('Unexpectedly large nonnumeric vector')
            if t==16:obj.value=[value(self.obj(path+'/string')) for j in range(n)]
            else:obj.value=[self.obj(path+f'/[{j}]') for j in range(n)]
        elif t==21:obj.value=self.bc(path)
        elif t==25:pass
        else:raise ValueError(f'Unsupported serialized type {t} at {self.pos} {path}')
        if flags&512:obj.attrs=self.obj(path+'/@')
        return obj
    def parse(self):
        magic=self.read(2)
        if magic!=b'X\n':raise ValueError(f'Only XDR supported: {magic}')
        self.header={'serialization_version':self.i(),'writer_version':self.i(),'minimum_reader_version':self.i()}
        if self.header['serialization_version']==3:self.header['encoding']=self.read(self.i()).decode()
        elif self.header['serialization_version']!=2:raise ValueError(self.header)
        root=self.obj('root')
        if self.f.read(1):raise ValueError('Unexpected trailing bytes')
        self.guard()
        return root

def main():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input');ap.add_argument('output');ap.add_argument('--keep-all-numeric',action='store_true')
    ap.add_argument('--rss-ceiling-gib',type=float,default=24.0);ap.add_argument('--disk-reserve-gib',type=float,default=75.0);a=ap.parse_args()
    p=Path(a.input);magic=p.open('rb').read(6)
    opener=gzip.open if magic[:2]==b'\x1f\x8b' else bz2.open if magic[:3]==b'BZh' else lzma.open if magic[:6]==b'\xfd7zXZ\x00' else open
    fixed={'rss_unit_multiplier':RSS_UNIT,'rss_ceiling_bytes':int(a.rss_ceiling_gib*1024**3),'disk_reserve_bytes':int(a.disk_reserve_gib*1024**3),'platform':platform.platform(),'archived_parser_sha256':ARCHIVED_SHA256,'compression':opener.__module__ if opener is not open else 'none'}
    with opener(p,'rb') as f:
        r=Reader(f,a.output,a.keep_all_numeric,fixed['rss_ceiling_bytes'],fixed['disk_reserve_bytes'])
        try:
            root=r.parse()
            with (r.out/'structure.pickle').open('wb') as s:pickle.dump(root,s,protocol=5)
            report={'status':'PARSED','header':r.header,'nodes':r.nodes,'uncompressed_bytes':r.pos,'elapsed_s':time.time()-r.start,'peak_RSS_bytes':rss_bytes(),**fixed,'root_attribute_names':list(attrs(root)),'payloads':r.payloads}
        except Exception as e:
            report={'status':'FAILED','error':repr(e),'nodes':r.nodes,'uncompressed_bytes':r.pos,'peak_RSS_bytes':rss_bytes(),**fixed,'payloads':r.payloads}
            raise
        finally:(r.out/'parse_report.json').write_text(json.dumps(report,indent=2))
    print({k:v for k,v in report.items() if k!='payloads'})
if __name__=='__main__':main()
