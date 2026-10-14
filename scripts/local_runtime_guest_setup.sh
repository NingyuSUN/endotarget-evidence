#!/usr/bin/env bash
set -euo pipefail
# Installs only into the newly imported project-owned distribution.
if [[ "${WSL_DISTRO_NAME:-}" != "EndoTargetRuntime" ]] || [[ "$(id -u)" != 0 ]]; then
  echo "Refuse setup outside root in EndoTargetRuntime" >&2
  exit 64
fi
project=/mnt/g/projects/te
runtime=/opt/endotarget-local
test -d "$project/.git"
mkdir -p "$runtime" /home/analyst/te
if [[ -e "$runtime/setup.complete" ]]; then
  echo "Preserve completed setup; inspect its receipt instead of reinstalling" >&2
  exit 65
fi
# Prevent package post-install scripts from launching a default networked daemon.
if [[ -e /usr/sbin/policy-rc.d ]]; then
  cp -p /usr/sbin/policy-rc.d "$runtime/policy-rc.d.original"
fi
printf '#!/bin/sh\nexit 101\n' > /usr/sbin/policy-rc.d
chmod 755 /usr/sbin/policy-rc.d
export DEBIAN_FRONTEND=noninteractive
apt-get update
packages=(docker.io python3 openjdk-21-jre-headless ca-certificates iptables procps xz-utils rsync git curl)
pins=()
: > "$runtime/apt-selected-versions.tsv"
for package in "${packages[@]}"; do
  version=$(apt-cache policy "$package" | awk '/Candidate:/ {print $2}')
  if [[ -z "$version" || "$version" == "(none)" ]]; then
    echo "Missing package candidate: $package" >&2
    exit 66
  fi
  printf '%s\t%s\n' "$package" "$version" >> "$runtime/apt-selected-versions.tsv"
  pins+=("$package=$version")
done
apt-get install -y --no-install-recommends "${pins[@]}"
dpkg-query -W -f='${Package}\t${Version}\n' > "$runtime/installed-packages.tsv"
# Match the prior task UID/GID inside this new project environment only.
if ! getent group analyst >/dev/null; then groupadd --gid 1004 analyst; fi
if ! id analyst >/dev/null 2>&1; then
  useradd --uid 1003 --gid 1004 --home-dir /home/analyst --shell /bin/bash analyst
fi
test "$(id -u analyst)" = 1003
test "$(id -g analyst)" = 1004
if ! mountpoint -q /home/analyst/te; then
  mount --bind "$project" /home/analyst/te
fi
test /home/analyst/te/README.md -ef "$project/README.md"
# No group access is granted to the privileged daemon; driver will be invoked explicitly.
if pgrep -x dockerd >/dev/null; then
  echo "Unexpected automatically started daemon; stop and inspect before continuing" >&2
  exit 67
fi
python3 - <<'PY'
import json,os,pathlib,subprocess,datetime
root=pathlib.Path('/home/analyst/te')
d={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
   'distribution':os.environ.get('WSL_DISTRO_NAME'),'execution_location':'LOCAL_PROJECT_WSL',
   'virtual_disk_host_location':r'G:\projects\te\work\local-runtime\distro',
   'project_mount':'/home/analyst/te','project_source':'/mnt/g/projects/te',
   'uid':1003,'gid':1004,'policy_rc_d_exit':101,'docker_daemon_started':False,
   'docker_client':subprocess.check_output(['docker','--version'],text=True).strip(),
   'java':subprocess.run(['java','-version'],capture_output=True,text=True).stderr.strip(),
   'python':subprocess.check_output(['python3','--version'],text=True).strip(),
   'existing_Ubuntu_or_Windows_security_modified':False,
   'container_execution':False}
(root/'receipts/local_runtime_guest_setup.json').write_text(json.dumps(d,indent=2)+'\n')
print(json.dumps(d,indent=2))
PY
touch "$runtime/setup.complete"
