#!/usr/bin/env bash
set -euo pipefail
if [[ "${WSL_DISTRO_NAME:-}" != EndoTargetRuntime ]] || [[ "$(id -u)" != 0 ]]; then
  echo "Refuse daemon launch outside root in EndoTargetRuntime" >&2; exit 64
fi
project=/mnt/g/projects/te
runtime=/opt/endotarget-local
test -f "$runtime/setup.complete"
test -d "$project/.git"
mkdir -p /home/analyst/te /var/lib/endotarget-docker "$runtime"
if ! mountpoint -q /home/analyst/te; then mount --bind "$project" /home/analyst/te; fi
test /home/analyst/te/README.md -ef "$project/README.md"
export DOCKER_HOST=unix:///run/endotarget-docker.sock
if docker info >/dev/null 2>&1 || pgrep -x dockerd >/dev/null; then
  echo "Existing Docker daemon detected; inspect rather than launch a duplicate" >&2; exit 65
fi
# No autostart: the setup policy also covers optional builder and inspection tools.
if ! command -v ip >/dev/null 2>&1; then
  version=$(apt-cache policy iproute2 | awk '/Candidate:/ {print $2}')
  test -n "$version"; test "$version" != '(none)'
  apt-get install -y --no-install-recommends "iproute2=$version"
  printf 'iproute2\t%s\n' "$version" > "$runtime/iproute2-selected-version.tsv"
fi
if ! docker buildx version >/dev/null 2>&1; then
  version=$(apt-cache policy docker-buildx | awk '/Candidate:/ {print $2}')
  test -n "$version"; test "$version" != '(none)'
  apt-get install -y --no-install-recommends "docker-buildx=$version"
  printf 'docker-buildx\t%s\n' "$version" > "$runtime/buildx-selected-version.tsv"
fi
python3 "$project/scripts/local_runtime_network_snapshot.py" "$runtime/network.before.json"
nohup dockerd \
  --host="$DOCKER_HOST" --group=root \
  --pidfile=/run/endotarget-docker.pid --exec-root=/run/endotarget-docker \
  --data-root=/var/lib/endotarget-docker --bridge=none \
  --iptables=false --ip6tables=false --ip-forward=false --ip-masq=false \
  --userland-proxy=false --exec-opt=native.cgroupdriver=cgroupfs \
  > "$project/logs/local_dockerd.log" 2>&1 &
daemon_pid=$!
for attempt in {1..60}; do
  if docker info >/dev/null 2>&1; then break; fi
  if ! kill -0 "$daemon_pid" 2>/dev/null; then tail -60 "$project/logs/local_dockerd.log"; exit 66; fi
  sleep 1
done
docker info >/dev/null
python3 "$project/scripts/local_runtime_network_snapshot.py" "$runtime/network.after.json"
python3 - <<'PY'
import json,pathlib,subprocess,datetime,os,hashlib
p=pathlib.Path('/mnt/g/projects/te'); r=pathlib.Path('/opt/endotarget-local')
before=json.loads((r/'network.before.json').read_text())
after=json.loads((r/'network.after.json').read_text())
changes={k:{'before':before[k],'after':after[k]} for k in before if before[k]!=after[k]}
info=json.loads(subprocess.check_output(['docker','info','--format','{{json .}}'],text=True))
pid=int(pathlib.Path('/run/endotarget-docker.pid').read_text())
args=pathlib.Path(f'/proc/{pid}/cmdline').read_bytes().decode().split('\0')
receipt={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'distribution':os.environ['WSL_DISTRO_NAME'],'daemon_pid':pid,'daemon_arguments':args,
 'docker_version':subprocess.check_output(['docker','version','--format','{{json .}}'],text=True).strip(),
 'buildx_version':subprocess.check_output(['docker','buildx','version'],text=True).strip(),
 'storage_driver':info['Driver'],'docker_root':info['DockerRootDir'],
 'security_options':info.get('SecurityOptions'),'cgroup_driver':info.get('CgroupDriver'),
 'socket':'unix:///run/endotarget-docker.sock','tcp_daemon_listener':False,
 'project_bind':'/home/analyst/te','project_source':'/mnt/g/projects/te',
 'network_before':before,'network_after':after,'network_changes':changes,
 'network_unchanged':not changes,'scientific_containers_executed':False,
 'script_sha256':{f:hashlib.sha256((p/f).read_bytes()).hexdigest() for f in
 ['scripts/local_runtime_start.sh','scripts/local_runtime_network_snapshot.py']}}
(p/'receipts/local_runtime_daemon.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'daemon_ready':True,'pid':pid,'network_unchanged':not changes,'changes':changes}))
if changes: raise SystemExit('Network state changed: inspect before proceeding')
PY
echo "Project daemon ready; this process keeps its WSL distribution active."
wait "$daemon_pid"
