#!/usr/bin/env python3
"""Build the additive local v3 recipe without changing the frozen 117-package lock."""
import datetime
import hashlib
import json
import os
import pathlib
import shutil
import subprocess

if os.environ.get('WSL_DISTRO_NAME') != 'EndoTargetRuntime' or os.geteuid() != 0:
    raise SystemExit('Refuse outside root in EndoTargetRuntime')
root = pathlib.Path('/home/analyst/te')
if not root.samefile('/mnt/g/projects/te'):
    raise SystemExit('Project bind missing')
env = dict(os.environ, DOCKER_HOST='unix:///run/endotarget-docker.sock')
subprocess.run(['docker', 'info'], env=env, check=True, stdout=subprocess.DEVNULL)
stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
context = pathlib.Path('/opt/endotarget-local/build-context-' + stamp)
context.mkdir()
sources = {'Dockerfile': 'containers/statistics-local/Dockerfile',
           '.dockerignore': 'containers/statistics-local/.dockerignore',
           'statistical-runtime.explicit.txt': 'requirements/statistical-runtime.explicit.txt',
           'os-packages.lock.tsv': 'containers/statistics-local/os-packages.lock.tsv'}
hashes = {}
for name, relative in sources.items():
    src, dest = root / relative, context / name
    shutil.copyfile(src, dest)
    hashes[relative] = hashlib.sha256(src.read_bytes()).hexdigest()
    assert hashlib.sha256(dest.read_bytes()).hexdigest() == hashes[relative]
tag = 'endotarget-statistics:local-v3-' + stamp.lower()
log = root / 'logs' / ('local_container_build_' + stamp + '.log')
cmd = ['docker', 'build', '--network=host', '--progress=plain', '--tag', tag, str(context)]
receipt = {'started_utc': stamp, 'distribution': 'EndoTargetRuntime', 'execution_location': 'LOCAL',
           'recipe_version': 'local_v3_add_procps_and_git', 'source_hashes': hashes, 'build_command': cmd,
           'build_context': str(context), 'image_tag': tag, 'build_log': str(log.relative_to(root)),
           'build_succeeded': False, 'nextflow_executed': False,
           'required_os_packages': dict(line.split('\t') for line in
             (root / sources['os-packages.lock.tsv']).read_text().splitlines())}
target = root / 'receipts/local_container_build.json'
if target.exists():
    raise SystemExit('Preserve previous build receipt before another build')
target.write_text(json.dumps(receipt, indent=2) + '\n')
with log.open('w') as stream:
    result = subprocess.run(cmd, env=env, stdout=stream, stderr=subprocess.STDOUT)
receipt['exit_code'] = result.returncode
receipt['build_log_sha256'] = hashlib.sha256(log.read_bytes()).hexdigest()
receipt['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
try:
    if result.returncode == 0:
        inspect = json.loads(subprocess.check_output(['docker', 'image', 'inspect', tag], env=env, text=True))[0]
        receipt.update(image_id=inspect['Id'], architecture=inspect['Architecture'], os=inspect['Os'],
                       size_bytes=inspect['Size'], image_config=inspect['Config'], repo_digests=inspect.get('RepoDigests', []))
        # Read the installed metadata with network disabled; no model execution.
        program = '''import glob,json,subprocess
packages=[]
for path in glob.glob('/opt/endo/conda-meta/*.json'):
 x=json.load(open(path));packages.append({k:x.get(k) for k in ('name','version','build','fn','md5','url')})
print(json.dumps({'scientific_packages':sorted(packages,key=lambda x:x['name']),
 'os_packages':subprocess.check_output(['dpkg-query','-W']+[line.split()[0] for line in open('/usr/local/share/endotarget-local-os-packages.lock.tsv')],text=True),
 'ps_version':subprocess.check_output(['ps','--version'],text=True),
 'git_version':subprocess.check_output(['git','--version'],text=True),
 'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd='/home/analyst/te',text=True).strip()}))
'''
        verify = ['docker', 'run', '--rm', '--network', 'none', '--memory', '2g', '--user', '1003:1004',
                  '--read-only', '--tmpfs', '/tmp:rw,nosuid,nodev,size=64m,mode=0700,uid=1003,gid=1004',
                  '-v', str(root)+':'+str(root)+':ro',
                  '--entrypoint', '/opt/endo/bin/python', tag, '-c', program]
        receipt['runtime_inspection_command'] = verify
        runtime = json.loads(subprocess.check_output(verify, env=env, text=True))
        receipt['runtime_inspection'] = runtime
        expected = {}
        for line in (root / sources['statistical-runtime.explicit.txt']).read_text().splitlines():
            if line.startswith('https://'):
                url, md5 = line.split('#')
                expected[url.rsplit('/', 1)[1]] = md5
        actual = {p['fn']: p['md5'] for p in runtime['scientific_packages']}
        if len(expected) != 117 or actual != expected:
            raise ValueError('Installed scientific package filenames/MD5 differ from frozen 117-package explicit lock')
        os_packages = {name.split(':', 1)[0]: version for name, version in
                       (line.split('\t') for line in runtime['os_packages'].splitlines())}
        if os_packages != receipt['required_os_packages']:
            raise ValueError('Pinned Debian metrics package versions differ')
        receipt.update(build_succeeded=True, runtime_inspection=runtime, scientific_lock_package_count=len(actual),
                       scientific_lock_exact_match=True, metrics_ps_available=True, provenance_git_available=True)
except Exception as error:
    receipt.update(build_succeeded=False, verification_error=repr(error))
receipt['script_sha256'] = hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
target.write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({k: receipt.get(k) for k in ('build_succeeded', 'exit_code', 'image_tag', 'image_id',
                                            'size_bytes', 'build_log', 'verification_error')}))
raise SystemExit(0 if receipt['build_succeeded'] else 1)
