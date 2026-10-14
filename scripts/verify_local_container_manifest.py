#!/usr/bin/env python3
"""Obtain and verify a real local registry manifest, then close the registry."""
import datetime,hashlib,json,os,pathlib,re,socket,subprocess,time,urllib.request
if os.environ.get("WSL_DISTRO_NAME")!="EndoTargetRuntime" or os.geteuid()!=0:
    raise SystemExit("Refuse outside root in EndoTargetRuntime")
root=pathlib.Path("/home/analyst/te")
if not root.samefile("/mnt/g/projects/te"): raise SystemExit("Project bind missing")
env=dict(os.environ,DOCKER_HOST="unix:///run/endotarget-docker.sock")
build=json.loads((root/"receipts/local_container_build.json").read_text())
if not build["build_succeeded"]: raise SystemExit("Build not successful")
stamp=datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
folder=root/"work/local-runtime"/("manifest-"+stamp); folder.mkdir()
receipt={"started_utc":stamp,"scope":"Temporary loopback registry; reachable from local host processes",
 "external_publication":False,"port":5500,"registry_stopped":False,"verified":False,"commands":[]}
target=root/"receipts/local_container_manifest.json"
if target.exists(): raise SystemExit("Preserve existing receipt before retry")
container=None
opener=urllib.request.build_opener(urllib.request.ProxyHandler({}))
def run(args):
    receipt["commands"].append(args)
    p=subprocess.run(args,env=env,capture_output=True,text=True)
    log=folder/(str(len(receipt["commands"]))+".log")
    log.write_text(p.stdout+p.stderr)
    if p.returncode: raise RuntimeError(f"{args[0:3]} exited {p.returncode}; see {log}")
    return p.stdout
def fetch(digest,label):
    url="http://127.0.0.1:5500/v2/endotarget-statistics/manifests/"+digest
    request=urllib.request.Request(url,headers={"Accept":", ".join([
        "application/vnd.oci.image.index.v1+json","application/vnd.oci.image.manifest.v1+json",
        "application/vnd.docker.distribution.manifest.v2+json",
        "application/vnd.docker.distribution.manifest.list.v2+json"])})
    with opener.open(request,timeout=30) as response:
        body=response.read(); header=response.headers.get("Docker-Content-Digest")
    calculated="sha256:"+hashlib.sha256(body).hexdigest()
    if calculated!=digest or header!=digest: raise ValueError("Manifest bytes/header/digest mismatch")
    file=folder/(label+".json"); file.write_bytes(body)
    return json.loads(body),{"requested_digest":digest,"header_digest":header,
         "raw_bytes_sha256":calculated,"file":str(file.relative_to(root))}
try:
    with socket.socket() as sock: sock.bind(("127.0.0.1",5500))
    run(["docker","pull","registry:3"])
    registry=json.loads(run(["docker","image","inspect","registry:3"]))[0]
    registry_ref=registry["RepoDigests"][0]
    receipt["registry_image_ref"]=registry_ref
    storage=pathlib.Path("/var/lib/endotarget-registry"); storage.mkdir(exist_ok=True)
    os.chown(storage,1003,1004)
    configuration=pathlib.Path("/opt/endotarget-local/registry.yml")
    configuration.write_text("version: 0.1\nlog:\n  level: info\nstorage:\n  filesystem:\n    rootdirectory: /var/lib/registry\nhttp:\n  addr: 127.0.0.1:5500\n")
    os.chmod(configuration,0o644)
    receipt["registry_configuration"] = configuration.read_text()
    name="endotarget-registry-"+stamp.lower()
    container=run(["docker","run","-d","--name",name,"--network=host",
      "--user","1003:1004","--read-only","--cap-drop=ALL","--security-opt=no-new-privileges",
      "--memory=256m","--cpus=1","--tmpfs","/tmp:rw,nosuid,nodev,size=64m,mode=0700,uid=1003,gid=1004",
      "-e","REGISTRY_HTTP_ADDR=127.0.0.1:5500","-e","OTEL_TRACES_EXPORTER=none",
      "-e","REGISTRY_HTTP_DEBUG_PROMETHEUS_ENABLED=false",
      "--mount","type=bind,source=/var/lib/endotarget-registry,target=/var/lib/registry",
      "--mount","type=bind,source=/opt/endotarget-local/registry.yml,target=/etc/distribution/local.yml,readonly",
      registry_ref,"/etc/distribution/local.yml"]).strip()
    receipt["registry_container"]=container
    for attempt in range(30):
        try:
            with opener.open("http://127.0.0.1:5500/v2/",timeout=2) as response:
                assert response.status==200
            break
        except Exception:
            if attempt==29: raise
            time.sleep(1)
    listeners=run(["ss","-H","-lntp"])
    relevant=[s for s in listeners.splitlines() if ":5500" in s or '"registry"' in s]
    if len(relevant)!=1 or "127.0.0.1:5500" not in relevant[0]:
        raise ValueError("Unexpected registry listener: "+repr(relevant))
    receipt["registry_listeners"]=relevant
    inspection=json.loads(run(["docker","inspect",container]))[0]
    receipt["registry_controls"]={k:inspection["HostConfig"][k] for k in
        ["NetworkMode","ReadonlyRootfs","CapDrop","SecurityOpt","Memory","Binds","Mounts"]}
    receipt["registry_user"]=inspection["Config"]["User"]
    local_tag="127.0.0.1:5500/endotarget-statistics:"+stamp.lower()
    recorded=json.loads(run(["docker","image","inspect",build["image_tag"]]))[0]
    if recorded["Id"]!=build["image_id"]: raise ValueError("Build tag no longer points to recorded image")
    run(["docker","tag",build["image_id"],local_tag])
    output=run(["docker","push",local_tag])
    digest=re.findall(r"digest: (sha256:[0-9a-f]{64})",output)[-1]
    parent,evidence=fetch(digest,"push-manifest")
    receipt["push_digest"]=digest; receipt["push_manifest"]=evidence
    if "manifests" in parent:
        candidates=[m for m in parent["manifests"] if m.get("platform",{}).get("os")=="linux"
                    and m.get("platform",{}).get("architecture")=="amd64"]
        if len(candidates)!=1: raise ValueError("Ambiguous amd64 manifest")
        digest=candidates[0]["digest"]
        manifest,evidence=fetch(digest,"platform-manifest")
        receipt["platform_manifest"]=evidence
    else: manifest=parent
    receipt["manifest_config_digest"]=manifest["config"]["digest"]
    # Check the image config bytes, not a locally invented IMAGE@config-digest reference.
    built=json.loads(run(["docker","image","inspect",build["image_id"]]))[0]
    config_digest=manifest["config"]["digest"]
    with opener.open("http://127.0.0.1:5500/v2/endotarget-statistics/blobs/"+config_digest,timeout=30) as response:
        config_bytes=response.read()
    if "sha256:"+hashlib.sha256(config_bytes).hexdigest()!=config_digest:
        raise ValueError("Image config blob digest mismatch")
    (folder/"image-config.json").write_bytes(config_bytes)
    config=json.loads(config_bytes)
    if config["architecture"]!="amd64" or config["os"]!="linux":
        raise ValueError("Image config platform differs")
    if config["rootfs"]["diff_ids"]!=built["RootFS"]["Layers"]:
        raise ValueError("Image layer configuration differs from built image")
    for key in ["User","Env","Cmd","Entrypoint","WorkingDir"]:
        if (config["config"].get(key) or None)!=(built["Config"].get(key) or None):
            raise ValueError("Image configuration differs: "+key)
    identities={receipt["push_digest"]:"pushed_manifest_or_index",digest:"platform_manifest",config_digest:"image_config"}
    if built["Id"] not in identities:
        raise ValueError("Built Docker ID does not belong to verified manifest/config chain")
    receipt["docker_id_kind"]=identities[built["Id"]]
    receipt["build_image_id"]=built["Id"]
    receipt["platform_digest"]=digest
    # Docker29 containerd image IDs can denote the OCI index rather than image config.
    # Preserve that verified parent reference; the selected amd64 child is checked above.
    reference="127.0.0.1:5500/endotarget-statistics@"+receipt["push_digest"]
    run(["docker","pull","--platform=linux/amd64",reference])
    pulled=json.loads(run(["docker","image","inspect",reference]))[0]
    if pulled["Id"]!=built["Id"] or pulled["Architecture"]!="amd64" or pulled["Os"]!="linux":
        raise ValueError("Pulled platform image differs from build")
    receipt.update(verified=True,image_ref=reference,image_id=pulled["Id"],image_config_digest=config_digest,
                   repo_digests=pulled["RepoDigests"])
except Exception as error:
    receipt["error"]=repr(error)
finally:
    if container:
        try:
            run(["docker","stop",container])
            state=json.loads(run(["docker","inspect",container]))[0]["State"]
            receipt["registry_stopped"]=not state["Running"]
        except Exception as error:
            receipt["cleanup_error"]=repr(error); receipt["verified"]=False
    receipt["finished_utc"]=datetime.datetime.now(datetime.timezone.utc).isoformat()
    receipt["artifacts"]={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in folder.iterdir() if p.is_file()}
    receipt["script_sha256"]=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
    target.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({k:receipt.get(k) for k in ["verified","image_ref","registry_stopped","error","cleanup_error"]}))
raise SystemExit(0 if receipt["verified"] and receipt["registry_stopped"] else 1)
