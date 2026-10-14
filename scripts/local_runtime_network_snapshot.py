#!/usr/bin/env python3
"""Record network-state digests without publishing local interface addresses."""
import hashlib,json,pathlib,subprocess,sys
def digest(value):
    return hashlib.sha256(value.encode()).hexdigest()
def command(args):
    p=subprocess.run(args,capture_output=True,text=True)
    # iptables-save emits timestamps in comments even when rules do not change.
    body="\n".join(s for s in p.stdout.splitlines() if not s.startswith("#"))
    return {"exit_code":p.returncode,"stdout_sha256":digest(body),"stderr_sha256":digest(p.stderr.strip())}
result={
 "ipv4_forward":pathlib.Path("/proc/sys/net/ipv4/ip_forward").read_text().strip(),
 "ipv6_forward":pathlib.Path("/proc/sys/net/ipv6/conf/all/forwarding").read_text().strip(),
 "routes":command(["ip","-j","route","show","table","all"]),
 "interfaces":command(["ip","-j","address","show"]),
 "iptables":command(["iptables-save"]),
 "ip6tables":command(["ip6tables-save"]),
}
pathlib.Path(sys.argv[1]).write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
