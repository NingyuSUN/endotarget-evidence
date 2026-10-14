#!/usr/bin/env python3
"""Run only preserved synthetic fixtures in local Docker; never run cohort models.

Execute as root in EndoTargetRuntime with an actual registry IMAGE@sha256 digest.
The receipt records engineering reproducibility, not new scientific validation.
"""
import argparse
import csv
import datetime
import hashlib
import json
import math
import os
import pathlib
import re
import subprocess
import threading
import time

ROOT = pathlib.Path('/home/analyst/te')
# Match the pinned upstream 25.10.2 launcher (lines 422-438): cache serialization
# needs these module accesses on Java21, even when the one.jar runs directly.
JAVA_MODULES = (
    'java.lang', 'java.io', 'java.nio', 'java.net', 'java.util',
    'java.util.concurrent.locks', 'java.util.concurrent.atomic', 'java.nio.file.spi',
    'sun.nio.ch', 'sun.nio.fs', 'sun.net.www.protocol.http', 'sun.net.www.protocol.https',
    'sun.net.www.protocol.ftp', 'sun.net.www.protocol.file', 'jdk.internal.misc',
    'jdk.internal.vm', 'java.util.regex',
)
NEXTFLOW_JVM_OPTIONS = ['-Xmx512m', '-Dfile.encoding=UTF-8'] + [
    '--add-opens=java.base/' + name + '=ALL-UNNAMED' for name in JAVA_MODULES
]
FIXTURES = {
    'discovery': ('20261013T100312Z-439d82d9', 'main.nf',
                  {'DONOR_DE', 'RANK_TARGETS', 'MATCHED_NULLS'}, 24),
    'replication': ('20261013T103604Z-replication-cb59eb61', 'replication.nf',
                    {'HECA_DONOR_DE', 'EVALUATE_REPLICATION'}, 10),
}


def sha(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def safe_path(root, name):
    path = pathlib.Path(name)
    if path.is_absolute() or '..' in path.parts:
        raise ValueError('Path must stay within project: ' + str(name))
    full = (root / path).resolve()
    if not full.is_relative_to(root.resolve()):
        raise ValueError('Path escapes project: ' + str(name))
    return full


def validate_image(image):
    if not re.fullmatch(r'[^\s@]+@sha256:[0-9a-f]{64}', image):
        raise ValueError('An immutable registry image digest IMAGE@sha256:DIGEST is required')


def check_image_binding(root, image, inspected, build):
    if image not in inspected.get('RepoDigests', []):
        raise ValueError('Image inspection does not confirm requested registry manifest digest')
    if build.get('build_succeeded') is not True or build.get('image_id') != inspected.get('Id'):
        raise ValueError('Registry image is not bound to the completed local build')
    if not build.get('source_hashes'):
        raise ValueError('Completed build has no source hashes')
    for name, digest in build['source_hashes'].items():
        if sha(safe_path(root, name)) != digest:
            raise ValueError('Container build source bytes changed: ' + name)


def compare_tsv(expected, actual, tolerance=1e-8):
    """Compare schema, row order, text and every finite number; no selected-row checks."""
    with expected.open(newline='') as a, actual.open(newline='') as b:
        left, right = list(csv.reader(a, delimiter='\t')), list(csv.reader(b, delimiter='\t'))
    if not left or not right or left[0] != right[0]:
        raise ValueError('TSV header mismatch: ' + str(actual))
    if len(left) != len(right):
        raise ValueError('TSV row count mismatch: ' + str(actual))
    numeric = 0
    for row_index, (old, new) in enumerate(zip(left[1:], right[1:]), 2):
        if len(old) != len(left[0]) or len(new) != len(left[0]):
            raise ValueError(f'TSV column count mismatch: {actual}:{row_index}')
        for col, (x, y) in enumerate(zip(old, new), 1):
            if x == y:
                continue
            try:
                xf, yf = float(x), float(y)
            except ValueError:
                raise ValueError(f'TSV cell mismatch: {actual}:{row_index}:{col}: {x!r} != {y!r}') from None
            if not (math.isfinite(xf) and math.isfinite(yf) and
                    math.isclose(xf, yf, rel_tol=tolerance, abs_tol=tolerance)):
                raise ValueError(f'TSV cell mismatch: {actual}:{row_index}:{col}: {x!r} != {y!r}')
            numeric += 1
    return {'rows': len(left) - 1, 'columns': len(left[0]),
            'numerically_equal_nonidentical_cells': numeric,
            'expected_sha256': sha(expected), 'actual_sha256': sha(actual)}


def compare_tree(expected, actual, expected_count):
    old = {p.relative_to(expected).as_posix(): p for p in expected.rglob('*.tsv')}
    new = {p.relative_to(actual).as_posix(): p for p in actual.rglob('*.tsv')}
    if len(old) != expected_count or set(old) != set(new):
        raise ValueError(f'TSV inventory mismatch: expected {expected_count}, native={sorted(old)}, actual={sorted(new)}')
    return {name: compare_tsv(old[name], new[name]) for name in sorted(old)}


def check_trace(path, stages, resume):
    with path.open(newline='') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    names = [r.get('name', '').split(' (', 1)[0] for r in rows]
    status = 'CACHED' if resume else 'COMPLETED'
    if len(rows) != len(stages) or set(names) != stages or len(set(names)) != len(names):
        raise ValueError('Trace does not contain exactly expected task stages')
    if any(r.get('status') != status or r.get('exit') != '0' for r in rows):
        raise ValueError('Trace requires every task ' + status + ' with exit 0')
    return rows


def clean_environment(source, run_dir, run_id):
    # Deliberate allowlist: no native conda, R/Python library, Java or Nextflow settings.
    env = {k: source[k] for k in ('WSL_DISTRO_NAME', 'WSL_INTEROP', 'WSLENV', 'LANG', 'LC_ALL', 'TERM') if k in source}
    env.update(PATH='/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin',
               HOME='/root', DOCKER_HOST='unix:///run/endotarget-docker.sock',
               ENDOTARGET_CONTAINER_UID='1003', ENDOTARGET_CONTAINER_GID='1004',
               ENDOTARGET_LOCAL_CHECK_ID=run_id, NXF_HOME=str(run_dir / 'nxf_home'),
               NXF_OFFLINE='true', NXF_DISABLE_CHECK_LATEST='true', NXF_ANSI_LOG='false',
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    return env


def check_container_inspect(item, root, workdir):
    cfg, host = item['Config'], item['HostConfig']
    if cfg.get('User') != '1003:1004':
        raise ValueError('Container did not use required non-root UID/GID')
    if host.get('NetworkMode') != 'none' or host.get('ReadonlyRootfs') is not True:
        raise ValueError('Container lacks network-none/read-only-root isolation')
    if host.get('Memory') != 2147483648 or host.get('Privileged') is not False:
        raise ValueError('Container memory or privilege settings differ')
    if host.get('SecurityOpt') and any('unconfined' in x for x in host['SecurityOpt']):
        raise ValueError('Container security profile was weakened')
    temporary = host.get('Tmpfs', {}).get('/tmp', '')
    if not all(x in temporary.split(',') for x in ('rw', 'nosuid', 'nodev', 'uid=1003', 'gid=1004')):
        raise ValueError('Container /tmp does not preserve the owned private tmpfs')
    mounts = {x['Destination']: x for x in item.get('Mounts', [])}
    if str(root) not in mounts or mounts[str(root)].get('RW') is not False:
        raise ValueError('Project is not mounted read-only')
    if str(workdir) not in mounts or mounts[str(workdir)].get('RW') is not True:
        raise ValueError('Task work directory is not mounted writable')
    if any('docker.sock' in x.get('Source', '') or 'docker.sock' in x['Destination'] for x in mounts.values()):
        raise ValueError('Container received a Docker daemon socket')


class ContainerObserver:
    """Inspect create events immediately because Nextflow removes task containers."""
    def __init__(self, env, directory, run_id):
        self.env, self.directory, self.items, self.errors = env, directory, [], []
        self.command = ['docker', 'events', '--filter', 'type=container', '--filter', 'event=create',
                        '--filter', 'label=endotarget.local-check=' + run_id, '--format', '{{json .}}']
        self.stderr = (directory / 'docker_events.stderr').open('w')
        self.process = subprocess.Popen(self.command, env=env, text=True, stdout=subprocess.PIPE, stderr=self.stderr)
        self.thread = threading.Thread(target=self.observe, daemon=True)
        self.thread.start()
        time.sleep(0.2)
        if self.process.poll() is not None:
            raise ValueError('Docker event observer failed to start')

    def observe(self):
        for line in self.process.stdout:
            try:
                event = json.loads(line)
                cid = event.get('id') or event['Actor']['ID']
                inspect = json.loads(subprocess.check_output(['docker', 'inspect', cid], env=self.env, text=True))[0]
                self.items.append({'event': event, 'inspect': inspect})
                (self.directory / ('container_' + cid + '.json')).write_text(json.dumps(self.items[-1], indent=2) + '\n')
            except Exception as error:
                self.errors.append(str(error))

    def close(self):
        self.process.terminate()
        self.process.wait(timeout=10)
        self.thread.join(timeout=10)
        self.stderr.close()
        if self.thread.is_alive():
            self.errors.append('Docker observer thread did not stop')


def snapshot_inputs(root):
    files = {root / n for n in ('main.nf', 'nextflow.config', 'replication.nf', 'replication.config',
                               'conf/statistics_container.config', 'conf/local_container.config',
                               'scripts/check_local_container_nextflow.py')}
    for fixture, _, _, _ in FIXTURES.values():
        folder = root / 'work/synthetic' / fixture
        files.add(folder / 'params.json')
        params = json.loads((folder / 'params.json').read_text())
        if params.get('synthetic_test') is not True:
            raise ValueError('Historical fixture must explicitly authorize only synthetic tests')
        manifest_name = params.get('evaluation_freeze') or params['pre_effect_freeze']
        manifest = safe_path(root, manifest_name)
        files.add(manifest)
        for name, digest in json.loads(manifest.read_text())['artifacts'].items():
            path = safe_path(root, name)
            if sha(path) != digest:
                raise ValueError('Frozen artifact changed before execution: ' + name)
            files.add(path)
        files.update(p for p in (folder / 'results_first').rglob('*') if p.is_file())
    return {str(p.relative_to(root)): sha(p) for p in sorted(files)}


def launch(root, run_dir, env, image, kind, resume):
    fixture, script, stages, count = FIXTURES[kind]
    folder = root / 'work/synthetic' / fixture
    label = kind + ('_resume' if resume else '_fresh')
    trace = run_dir / (label + '_trace.tsv')
    out = run_dir / (label + '_results')
    cache = run_dir / (kind + '_work')
    if out.exists():
        raise ValueError('Refuse to overwrite output: ' + str(out))
    configs = ['replication.config'] if kind == 'replication' else ['nextflow.config']
    configs += ['conf/statistics_container.config', 'conf/local_container.config']
    command = ['java'] + NEXTFLOW_JVM_OPTIONS + ['-jar', str(root / 'work/local-runtime/nextflow/nextflow-25.10.2-one.jar'),
               '-C', ','.join(str(root / name) for name in configs),
               '-log', str(run_dir / (label + '_nextflow.log')), 'run', str(root / script),
               '-profile', 'docker,statistics_container,local_container',
               '-params-file', str(folder / 'params.json'), '--runtime_mode', 'docker',
               '--container_image', image, '--outdir', str(out),
               '-work-dir', str(cache), '-with-trace', str(trace), '-dump-hashes', 'json']
    name = 'local_' + run_dir.name + '_' + kind
    command += ['-resume', name] if resume else ['-name', name]
    log = run_dir / (label + '_console.log')
    (run_dir / (label + '_command.json')).write_text(json.dumps(command, indent=2) + '\n')
    with log.open('w') as stream:
        result = subprocess.run(command, cwd=run_dir, env=env, stdout=stream, stderr=subprocess.STDOUT, timeout=900)
    if result.returncode:
        raise ValueError(f'{label} Nextflow exit {result.returncode}; see {log}')
    rows = check_trace(trace, stages, resume)
    comparison = compare_tree(folder / 'results_first', out, count)
    if resume:
        for path in (run_dir / (kind + '_fresh_results')).rglob('*'):
            if path.is_file() and sha(path) != sha(out / path.relative_to(run_dir / (kind + '_fresh_results'))):
                raise ValueError('Cached output bytes changed: ' + str(path))
    tasks = []
    for row in rows:
        work = pathlib.Path(row['workdir']).resolve()
        if not work.is_relative_to(cache.resolve()):
            raise ValueError('Trace task work directory escapes isolated cache')
        run = work / '.command.run'
        shell = work / '.command.sh'
        text = run.read_text()
        if 'docker run' not in text or image not in text:
            raise ValueError('Task launcher does not invoke the immutable Docker image')
        if 'source ' in shell.read_text() and 'statistical_runtime_env.sh' in shell.read_text():
            raise ValueError('Task sourced the historical native runtime')
        tasks.append({'name': row['name'], 'workdir': str(work),
                      'command_run': str(run.relative_to(root)), 'command_run_sha256': sha(run),
                      'command_sh': str(shell.relative_to(root)), 'command_sh_sha256': sha(shell)})
    return {'pipeline': kind, 'resume': resume, 'command': command, 'exit_code': result.returncode,
            'console': str(log.relative_to(root)), 'trace': str(trace.relative_to(root)),
            'tasks': tasks, 'trace_rows': rows, 'all_scientific_tsv_comparisons': comparison,
            'output': str(out.relative_to(root)), 'cached_output_bytes_identical': True if resume else None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--image', required=True, help='Actual registry IMAGE@sha256:DIGEST, never a local image ID')
    args = parser.parse_args()
    validate_image(args.image)
    if os.environ.get('WSL_DISTRO_NAME') != 'EndoTargetRuntime' or os.geteuid() != 0:
        raise SystemExit('Refuse outside root in EndoTargetRuntime')
    if not ROOT.samefile('/mnt/g/projects/te') or pathlib.Path(__file__).resolve().parent.parent != ROOT.resolve():
        raise SystemExit('Expected local project bind /mnt/g/projects/te -> /home/analyst/te')
    receipt_path = ROOT / 'receipts/local_container_nextflow.json'
    if receipt_path.exists():
        raise SystemExit('Preserve existing local container receipt; choose review before rerun')
    original = snapshot_inputs(ROOT)
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    run_id = stamp + '-' + os.urandom(4).hex()
    # Frozen synthetic evaluation enforces this location for external model outputs.
    run_dir = ROOT / 'work/synthetic/local-container' / run_id
    run_dir.mkdir(parents=True)
    env = clean_environment(os.environ, run_dir, run_id)
    image_info = json.loads(subprocess.check_output(['docker', 'image', 'inspect', args.image], env=env, text=True))[0]
    build_path = ROOT / 'receipts/local_container_build.json'
    build = json.loads(build_path.read_text())
    check_image_binding(ROOT, args.image, image_info, build)
    java = subprocess.run(['java', '-version'], env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True).stdout
    if not re.search(r'version "21[.\"]', java):
        raise SystemExit('Host Java 21 is required')
    receipt = {'status': 'RUNNING', 'execution_location': 'LOCAL_EndoTargetRuntime',
               'data_kind': 'synthetic_test_only', 'scientific_validation': False,
               'real_cohort_effects_computed': False, 'AWS_execution_verified': False,
               'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'image': args.image, 'image_id': image_info['Id'], 'image_inspect': image_info,
               'build_receipt': str(build_path.relative_to(ROOT)), 'build_receipt_sha256': sha(build_path),
               'build_image_id': build['image_id'],
               'host_java': java, 'java_heap': '512m', 'task_concurrency': 1,
               'java_launcher_source': 'https://raw.githubusercontent.com/nextflow-io/nextflow/v25.10.2/nextflow',
               'java_launcher_flags': NEXTFLOW_JVM_OPTIONS,
               'comparison_tolerance': {'relative': 1e-8, 'absolute': 1e-8},
               'run_directory': str(run_dir.relative_to(ROOT)), 'input_sha256_before': original,
               'environment': env, 'launches': []}
    observer = None
    try:
        observer = ContainerObserver(env, run_dir, run_id)
        for kind in FIXTURES:
            for resume in (False, True):
                receipt['launches'].append(launch(ROOT, run_dir, env, args.image, kind, resume))
                print(kind + (' resume CACHED' if resume else ' fresh COMPLETED'), flush=True)
        observer.close()
        receipt['docker_event_command'] = observer.command
        receipt['container_inspections'] = observer.items
        if observer.errors or len(observer.items) != 5:
            raise ValueError(f'Expected 5 actual containers, got {len(observer.items)}; observer errors={observer.errors}')
        taskdirs = {task['workdir'] for item in receipt['launches'] if not item['resume'] for task in item['tasks']}
        observed = set()
        for item in observer.items:
            inspect = item['inspect']
            if inspect['Config']['Image'] != args.image or inspect['Image'] != image_info['Id']:
                raise ValueError('Actual task container image differs from inspected digest')
            matches = taskdirs.intersection(m['Destination'] for m in inspect.get('Mounts', []))
            if len(matches) != 1:
                raise ValueError('Cannot map actual container to exactly one task directory')
            work = next(iter(matches))
            check_container_inspect(inspect, ROOT, pathlib.Path(work))
            observed.add(work)
        if observed != taskdirs:
            raise ValueError('Actual container evidence does not cover all 5 fresh tasks')
        receipt.update(status='PASS', docker_execution_verified=True,
                       fresh_tasks_completed=5, resumed_tasks_cached=5,
                       discovery_tsvs_compared=24, replication_tsvs_compared=10,
                       original_inputs_and_native_outputs_preserved=True)
    except Exception as error:
        receipt.update(status='FAIL', docker_execution_verified=False, error=str(error))
    finally:
        if observer is not None and observer.process.poll() is None:
            observer.close()
        # Hash the original inventory directly: a changed freeze must still yield a failure receipt.
        receipt['input_sha256_after'] = {}
        for name in original:
            try:
                receipt['input_sha256_after'][name] = sha(ROOT / name)
            except OSError as error:
                receipt['input_sha256_after'][name] = 'UNREADABLE: ' + str(error)
        if observer is not None:
            receipt['docker_event_command'] = observer.command
            receipt['container_inspections'] = observer.items
            receipt['container_observer_errors'] = observer.errors
        if receipt['input_sha256_after'] != original:
            receipt.update(status='FAIL', docker_execution_verified=False,
                           original_inputs_and_native_outputs_preserved=False, error='Input/native evidence bytes changed')
        receipt['completed_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        with receipt_path.open('x') as stream:
            json.dump(receipt, stream, indent=2, allow_nan=False)
            stream.write('\n')
    print(json.dumps({'status': receipt['status'], 'receipt': str(receipt_path), 'error': receipt.get('error')}))
    return 0 if receipt['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
