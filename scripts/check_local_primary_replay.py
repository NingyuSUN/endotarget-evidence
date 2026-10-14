#!/usr/bin/env python3
"""Bounded clean-container replay of the two frozen real primary analyses.

This executes five fresh tasks once. It does not import single-cell data, rerun
bootstrap/LODO, select new external lists, or establish new scientific validity.
"""
import argparse
import csv
import datetime
import itertools
import json
import math
import os
import pathlib
import re
import subprocess

from check_local_container_nextflow import (
    ContainerObserver, NEXTFLOW_JVM_OPTIONS, check_container_inspect, check_image_binding, check_trace,
    clean_environment, sha, validate_image,
)

ROOT = pathlib.Path('/home/analyst/te')
TOLERANCE = 1e-8  # Predeclared absolute tolerance only; never magnitude-scaled.
FREEZES = {'statistical_pre_effect_freeze': (43, 'STATISTICAL_IMPLEMENTATION_FROZEN_PRE_EFFECT'),
           'discovery_rank_freeze': (74, 'DISCOVERY_FROZEN_PRE_EXTERNAL'),
           'replication_implementation_freeze': (95, 'REPLICATION_IMPLEMENTATION_FROZEN_PRE_EXTERNAL')}
DE_FILES = {'de/' + n for n in ('de_all.tsv', 'de_fitted.tsv', 'design.tsv',
            'expression_filter.tsv', 'normalization.tsv', 'model.rds', 'provenance.json', 'sessionInfo.txt')}
DISCOVERY_FILES = DE_FILES | {'ranking/' + n for n in (
    'B0.tsv', 'B1.tsv', 'G1.tsv', 'B0_top20.tsv', 'B1_top20.tsv', 'G1_top20.tsv',
    'B0_top50.tsv', 'B1_top50.tsv', 'G1_top50.tsv', 'excluded_genes.tsv',
    'induced_edges.tsv', 'universe.tsv', 'provenance.json')} | {'nulls/' + n for n in (
    'B0_null_membership.tsv', 'B1_null_membership.tsv', 'G1_null_membership.tsv',
    'all_U_directions.tsv', 'bin_demands.tsv', 'boundaries.json', 'null_memberships.tsv',
    'status.json', 'universe_bins.tsv')}
HECA_FILES = DE_FILES | {'evaluation/' + n for n in (
    'bootstrap_interface.tsv', 'dossier_pool.tsv', 'dossier_selection.tsv', 'null_sets.tsv',
    'null_summary.json', 'provenance.json', 'selected_genes.tsv', 'summary.json')}
# Unknown fields default to exact comparison, so numeric IDs and new count fields
# cannot silently acquire floating tolerance. Discrete design indicators stay exact.
CONTINUOUS_COLUMNS = frozenset(('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val',
    'moderated_SE', 'df_total', 'CI_low', 'CI_high', 'null_matching_logCPM',
    'age_centered_years', 'norm_factor',
    'effective_library_size', 'score', 'weighted_degree', 'weight', 'external_logFC',
    'external_t', 'external_p', 'external_q', 'external_CI_low', 'external_CI_high',
    'signed_external_t', 'concordance', 'coverage', 'missing_lower', 'missing_upper'))
JOBS = {
    'discovery': {'script': 'main.nf', 'config': 'nextflow.config', 'stages':
        {'DONOR_DE', 'RANK_TARGETS', 'MATCHED_NULLS'}, 'source': 'work/statistical/discovery_primary',
        'files': DISCOVERY_FILES, 'tested_genes': 17733},
    'heca': {'script': 'replication.nf', 'config': 'replication.config', 'stages':
        {'HECA_DONOR_DE', 'EVALUATE_REPLICATION'}, 'source': 'work/statistical/heca_primary',
        'files': HECA_FILES, 'tested_genes': 17643},
}


def project_path(root, name):
    p = pathlib.Path(name)
    if p.is_absolute() or '..' in p.parts or p.as_posix() != name:
        raise ValueError('Require normalized project-relative path: ' + str(name))
    full = root / p
    if not full.resolve().is_relative_to(root.resolve()):
        raise ValueError('Path escapes project: ' + name)
    return full


def inventory(folder, expected):
    found = {p.relative_to(folder).as_posix(): p for p in folder.rglob('*') if p.is_file()}
    if set(found) != expected or any(p.is_symlink() for p in folder.rglob('*')):
        raise ValueError(f'Output inventory mismatch/symlink: {folder}; '
                         f'missing={sorted(expected-set(found))}, extra={sorted(set(found)-expected)}')
    return found


def compare_tsv(expected, actual):
    adjusted, rows, largest = 0, 0, 0.0
    with expected.open(newline='') as a, actual.open(newline='') as b:
        old, new = csv.reader(a, delimiter='\t'), csv.reader(b, delimiter='\t')
        header = next(old, None)
        if not header or header != next(new, None) or len(set(header)) != len(header):
            raise ValueError('TSV header mismatch: ' + str(actual))
        for rows, pair in enumerate(itertools.zip_longest(old, new), 1):
            left, right = pair
            if left is None or right is None or len(left) != len(header) or len(right) != len(header):
                raise ValueError(f'TSV row/column count mismatch: {actual}:{rows+1}')
            for field, x, y in zip(header, left, right):
                if x == y:
                    continue
                if field in CONTINUOUS_COLUMNS:
                    try:
                        xf, yf = float(x), float(y)
                    except ValueError:
                        xf = yf = math.nan
                    if math.isfinite(xf) and math.isfinite(yf) and math.isclose(xf, yf, rel_tol=0, abs_tol=TOLERANCE):
                        adjusted += 1
                        largest = max(largest, abs(xf-yf))
                        continue
                raise ValueError(f'TSV mismatch: {actual}:{rows+1}:{field}: {x!r} != {y!r}')
    return {'rows': rows, 'columns': len(header), 'numerically_equal_nonidentical_cells': adjusted,
            'maximum_absolute_difference': largest, 'expected_sha256': sha(expected), 'actual_sha256': sha(actual)}


def compare_json_values(old, new, pointer='$'):
    if isinstance(old, dict):
        if not isinstance(new, dict) or set(old) != set(new):
            raise ValueError('JSON schema differs: ' + pointer)
        for key in old:
            compare_json_values(old[key], new[key], pointer + '.' + key)
    elif isinstance(old, list):
        if not isinstance(new, list) or len(old) != len(new):
            raise ValueError('JSON list differs: ' + pointer)
        for i, (x, y) in enumerate(zip(old, new)):
            compare_json_values(x, y, pointer + f'[{i}]')
    elif isinstance(old, float):
        if isinstance(new, bool) or not isinstance(new, (int, float)) or not (
                math.isfinite(old) and math.isfinite(new) and math.isclose(old, new, rel_tol=0, abs_tol=TOLERANCE)):
            raise ValueError('JSON continuous value differs: ' + pointer)
    elif type(old) is not type(new) or old != new:
        raise ValueError('JSON discrete/text value differs: ' + pointer)


def verify_output_hashes(value, folder):
    """Validate generated metadata against its own files before comparing semantics."""
    for name, meta in value.items():
        digest = meta.get('sha256') if isinstance(meta, dict) else meta
        path = folder / pathlib.Path(name).name
        if not re.fullmatch('[0-9a-f]{64}', str(digest)) or sha(path) != digest:
            raise ValueError('Output provenance hash mismatch: ' + str(path))
        if isinstance(meta, dict) and 'bytes' in meta and path.stat().st_size != meta['bytes']:
            raise ValueError('Output provenance byte count mismatch: ' + str(path))


def generated_input_path(output, name):
    recorded = pathlib.Path(name)
    if not recorded.is_absolute() or not recorded.is_relative_to(ROOT):
        raise ValueError('Generated provenance input escapes canonical project: ' + name)
    roots = [p for p in (output, *output.parents) if (p / 'main.nf').is_file()]
    if not roots:
        raise ValueError('Cannot locate project for generated provenance input')
    return project_path(roots[0], recorded.relative_to(ROOT).as_posix())


def normalized_json(path, output):
    """Exclude only declared execution metadata; retain scientific parameters."""
    value = json.loads(path.read_text())
    name = path.relative_to(output).as_posix()
    if name == 'de/provenance.json':
        verify_output_hashes(value.pop('output_sha256'), output / 'de')
        for key in ('started_utc', 'completed_utc'):
            value.pop(key)
    elif name in ('ranking/provenance.json', 'nulls/status.json'):
        verify_output_hashes(value.pop('outputs'), path.parent)
        value.pop('finished_utc' if name.startswith('ranking') else 'completed_utc')
        # Verify every referenced input. Source inputs keep path/hash exact; only
        # task-local generated paths/hashes map to previously compared outputs.
        for role, item in value['inputs'].items():
            if name.startswith('ranking') and role in ('discovery_de', 'discovery_de_all'):
                relative = 'de/' + ('de_fitted.tsv' if role == 'discovery_de' else 'de_all.tsv')
            elif name.startswith('nulls') and role == 'discovery_de':
                relative = 'de/de_fitted.tsv'
            elif name.startswith('nulls') and role not in ('configuration', 'implementation'):
                relative = 'ranking/' + role
            else:
                continue
            reference = output / relative
            if (sha(reference) != item['sha256'] or pathlib.Path(item['path']).name != reference.name
                    or sha(generated_input_path(output, item['path'])) != item['sha256']):
                raise ValueError('Generated input provenance differs: ' + role)
            item.update(path=relative, sha256='SEMANTICALLY_COMPARED_OUTPUT')
        value['software'].pop('executable')
        if name.startswith('ranking'):
            for key, relative in [('universe_sha256', 'ranking/universe.tsv'), ('induced_graph_sha256', 'ranking/induced_edges.tsv')]:
                if value[key] != sha(output / relative):
                    raise ValueError('Ranking hash reference differs: ' + key)
                value[key] = 'SEMANTICALLY_COMPARED_OUTPUT'
    elif name == 'evaluation/provenance.json':
        verify_output_hashes(value.pop('output_sha256'), path.parent)
        verify_output_hashes(value['validated_external'].pop('outputs'), output / 'de')
        if value.pop('external_provenance_sha256') != sha(output / 'de/provenance.json'):
            raise ValueError('Evaluation did not bind actual external DE provenance')
        # Code hashes and freeze hashes remain exact. Git HEAD and wall time are
        # host execution metadata; source file hashes establish implementation.
        value.pop('completed_utc')
        value.pop('git_head')
    return value


def compare_outputs(expected, actual, job):
    old = inventory(expected, job['files'])
    new = inventory(actual, job['files'])
    tables = {n: compare_tsv(old[n], new[n]) for n in sorted(old) if n.endswith('.tsv')}
    semantic = {}
    for n in sorted(old):
        if n.endswith('.json'):
            compare_json_values(normalized_json(old[n], expected), normalized_json(new[n], actual), n)
            semantic[n] = {'scientific_semantics_equal': True, 'expected_sha256': sha(old[n]), 'actual_sha256': sha(new[n])}
    provenance = json.loads(new['de/provenance.json'].read_text())
    if any(provenance.get(k) != v for k, v in {'tested_family_size': job['tested_genes'],
            'estimable_genes': job['tested_genes'], 'significant_BH05': 0, 'data_kind': 'public_cohort',
            'scientific_effect_run': True}.items()):
        raise ValueError('Real primary DE expected cohort/gene-count/BH invariant failed')
    if tables['de/de_fitted.tsv']['rows'] != job['tested_genes']:
        raise ValueError('Real primary DE tested gene count failed')
    runtime = {}
    for n in ('de/model.rds', 'de/sessionInfo.txt'):
        if new[n].stat().st_size == 0:
            raise ValueError('Missing runtime artifact: ' + n)
        runtime[n] = {'expected_sha256': sha(old[n]), 'actual_sha256': sha(new[n]),
                      'comparison': 'record_only_native_serialized_model_or_host_session_metadata',
                      'bytes_equal': sha(old[n]) == sha(new[n])}
    return {'all_scientific_tsv_comparisons': tables, 'scientific_json_comparisons': semantic,
            'runtime_artifacts': runtime, 'all_output_files': sorted(new),
            'tested_genes': job['tested_genes'], 'BH_q_lt_05': 0}


def snapshot_inputs(root):
    files = set()
    frozen = {}
    for name, (count, phase) in FREEZES.items():
        p = root / 'receipts' / (name + '.json')
        value = json.loads(p.read_text())
        if value.get('phase') != phase or len(value.get('artifacts', {})) != count:
            raise ValueError('Freeze count/phase changed: ' + name)
        files.add(p)
        frozen[name] = {'artifact_count': count, 'sha256': sha(p)}
        for name, digest in value['artifacts'].items():
            path = project_path(root, name)
            if sha(path) != digest:
                raise ValueError('Frozen artifact changed before replay: ' + name)
            files.add(path)
    for job in JOBS.values():
        files.update(inventory(root / job['source'], job['files']).values())
    for name in ('scripts/check_local_primary_replay.py', 'scripts/check_local_container_nextflow.py',
                 'conf/local_container.config', 'conf/statistics_container.config',
                 'config/replication_run.json', 'receipts/local_container_build.json',
                 'receipts/local_container_nextflow.json', 'work/local-runtime/nextflow/nextflow-25.10.2-one.jar'):
        files.add(project_path(root, name))
    # Bind each historical DE artifact to recorded output hashes, including files
    # not themselves listed in the 74/95-entry freezes (e.g. HECA native results).
    for job in JOBS.values():
        source = root / job['source']
        for name in job['files']:
            if name.endswith('.json'):
                normalized_json(source / name, source)
    heca = json.loads((root / 'receipts/heca_primary_execution.json').read_text())
    files.add(root / 'receipts/heca_primary_execution.json')
    for name, digest in heca['output_sha256'].items():
        if sha(project_path(root, name)) != digest:
            raise ValueError('Native HECA receipt artifact differs: ' + name)
    return {p.relative_to(root).as_posix(): sha(p) for p in sorted(files)}, frozen


def launch(root, run_dir, env, image, kind):
    job = JOBS[kind]
    out, cache = run_dir / (kind + '_results'), run_dir / (kind + '_work')
    if out.exists() or cache.exists():
        raise ValueError('Refuse existing replay output/cache')
    # Discovery has no run params file in the historical run. Its frozen config
    # supplies all inputs. HECA preserves original discovery_freeze, never replay.
    params = ({'pre_effect_freeze': 'receipts/statistical_pre_effect_freeze.json'} if kind == 'discovery'
              else json.loads((root / 'config/replication_run.json').read_text()))
    params.update(runtime_mode='docker', container_image=image, outdir=str(out))
    if params.get('synthetic_test', False) is not False:
        raise ValueError('Real replay cannot use synthetic fixtures')
    params_file = run_dir / (kind + '_params.json')
    params_file.write_text(json.dumps(params, indent=2) + '\n')
    configs = [job['config'], 'conf/statistics_container.config', 'conf/local_container.config']
    trace = run_dir / (kind + '_trace.tsv')
    command = ['java'] + NEXTFLOW_JVM_OPTIONS + ['-jar', str(root / 'work/local-runtime/nextflow/nextflow-25.10.2-one.jar'),
               '-C', ','.join(str(root / n) for n in configs), '-log', str(run_dir / (kind + '_nextflow.log')),
               'run', str(root / job['script']), '-profile', 'docker,statistics_container,local_container',
               '-params-file', str(params_file), '-work-dir', str(cache), '-with-trace', str(trace),
               '-name', 'primary_' + run_dir.name + '_' + kind]
    console = run_dir / (kind + '_console.log')
    (run_dir / (kind + '_command.json')).write_text(json.dumps(command, indent=2) + '\n')
    with console.open('x') as stream:
        result = subprocess.run(command, cwd=run_dir, env=env, stdout=stream, stderr=subprocess.STDOUT, timeout=900)
    if result.returncode:
        raise ValueError(f'{kind} Nextflow exit {result.returncode}; see {console}')
    rows = check_trace(trace, job['stages'], False)
    tasks = []
    for row in rows:
        work = pathlib.Path(row['workdir']).resolve()
        if not work.is_relative_to(cache.resolve()):
            raise ValueError('Task escaped fresh replay work directory')
        launcher, shell = work / '.command.run', work / '.command.sh'
        if 'docker run' not in launcher.read_text() or image not in launcher.read_text():
            raise ValueError('Actual task launcher lacks pinned Docker image')
        if 'statistical_runtime_env.sh' in shell.read_text():
            raise ValueError('Actual task reused native runtime')
        tasks.append({'name': row['name'], 'workdir': str(work),
                      'command_run_sha256': sha(launcher), 'command_sh_sha256': sha(shell)})
    comparison = compare_outputs(root / job['source'], out, job)
    return {'pipeline': kind, 'resume': False, 'command': command, 'effective_params': params,
            'source': job['source'], 'output': str(out.relative_to(root)), 'trace_rows': rows,
            'tasks': tasks, 'exit_code': result.returncode, **comparison}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--image', required=True)
    args = parser.parse_args()
    validate_image(args.image)
    if os.environ.get('WSL_DISTRO_NAME') != 'EndoTargetRuntime' or os.geteuid() != 0:
        raise SystemExit('Refuse outside root in EndoTargetRuntime')
    if not ROOT.samefile('/mnt/g/projects/te') or pathlib.Path(__file__).resolve().parent.parent != ROOT.resolve():
        raise SystemExit('Expected EndoTargetRuntime local G: project bind')
    receipt_path = ROOT / 'receipts/local_primary_container_replay.json'
    if receipt_path.exists():
        raise SystemExit('Preserve existing real replay receipt; no automatic repeat')
    synthetic_path = ROOT / 'receipts/local_container_nextflow.json'
    synthetic = json.loads(synthetic_path.read_text())
    if synthetic.get('status') != 'PASS' or synthetic.get('docker_execution_verified') is not True or synthetic.get('image') != args.image:
        raise SystemExit('Actual same-image synthetic Docker gate must PASS before real replay')
    original, freezes = snapshot_inputs(ROOT)
    run_id = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + os.urandom(4).hex()
    run_dir = ROOT / 'work/local-primary-container' / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    env = clean_environment(os.environ, run_dir, run_id)
    image_info = json.loads(subprocess.check_output(['docker', 'image', 'inspect', args.image], env=env, text=True))[0]
    build_path = ROOT / 'receipts/local_container_build.json'
    build = json.loads(build_path.read_text())
    check_image_binding(ROOT, args.image, image_info, build)
    java = subprocess.run(['java', '-version'], env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True).stdout
    if not re.search(r'version "21[.\"]', java):
        raise SystemExit('Host Java 21 required')
    receipt = {'status': 'RUNNING', 'execution_location': 'LOCAL_EndoTargetRuntime', 'data_kind': 'public_cohort',
               'scope': 'bounded_frozen_real_primary_replay_3_discovery_plus_2_HECA_tasks',
               'new_scientific_validation': False, 'AWS_execution_verified': False,
               'bootstrap_500_rerun': False, 'LODO_22_rerun': False, 'raw_import_rerun': False,
               'HECA_selection_source': 'original_frozen_discovery_not_new_replay',
               'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'image': args.image, 'image_id': image_info['Id'], 'image_inspect': image_info,
               'build_receipt_sha256': sha(build_path), 'synthetic_receipt_sha256': sha(synthetic_path),
               'comparison_tolerance': {'relative': 0, 'absolute': TOLERANCE, 'continuous_columns': sorted(CONTINUOUS_COLUMNS),
                   'discrete_ids_counts_ranks_directions_filters': 'exact'},
               'runtime_artifact_exceptions': ['RDS serialization and sessionInfo host descriptions recorded, not byte-equality required',
                   'JSON timestamps/work paths/software executable/Git HEAD recorded; scientific values retained and referenced hashes checked'],
               'run_directory': str(run_dir.relative_to(ROOT)), 'input_sha256_before': original,
               'freezes': freezes, 'environment': env, 'host_java': java, 'java_heap': '512m',
               'task_concurrency': 1, 'launches': []}
    observer = None
    try:
        observer = ContainerObserver(env, run_dir, run_id)
        for kind in JOBS:
            receipt['launches'].append(launch(ROOT, run_dir, env, args.image, kind))
            print(kind + ' fresh real primary COMPLETED and compared', flush=True)
        observer.close()
        if observer.errors or len(observer.items) != 5:
            raise ValueError(f'Need exactly 5 actual containers: {len(observer.items)}; {observer.errors}')
        taskdirs = {task['workdir'] for run in receipt['launches'] for task in run['tasks']}
        observed = []
        for item in observer.items:
            value = item['inspect']
            if value['Config']['Image'] != args.image or value['Image'] != image_info['Id']:
                raise ValueError('Actual task image differs from completed build')
            matches = taskdirs.intersection(m['Destination'] for m in value.get('Mounts', []))
            if len(matches) != 1:
                raise ValueError('Cannot map container to exactly one fresh task')
            work = next(iter(matches))
            check_container_inspect(value, ROOT, pathlib.Path(work))
            observed.append(work)
        if len(set(observed)) != 5 or set(observed) != taskdirs:
            raise ValueError('Docker inspections do not cover exactly five fresh real tasks')
        receipt.update(status='PASS', docker_execution_verified=True, fresh_tasks_completed=5,
                       discovery_tsvs_compared=24, HECA_tsvs_compared=10,
                       real_primary_outputs_reproduced=True, original_inputs_and_native_outputs_preserved=True)
    except Exception as error:
        receipt.update(status='FAIL', docker_execution_verified=False, real_primary_outputs_reproduced=False, error=str(error))
    finally:
        if observer is not None:
            if observer.process.poll() is None:
                observer.close()
            receipt.update(docker_event_command=observer.command, container_inspections=observer.items,
                           container_observer_errors=observer.errors)
        after = {}
        for name in original:
            try:
                after[name] = sha(ROOT / name)
            except OSError as error:
                after[name] = 'UNREADABLE: ' + str(error)
        receipt['input_sha256_after'] = after
        if after != original:
            receipt.update(status='FAIL', real_primary_outputs_reproduced=False,
                           original_inputs_and_native_outputs_preserved=False, error='Protected source/native evidence changed')
        receipt['completed_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        with receipt_path.open('x') as stream:
            json.dump(receipt, stream, indent=2, allow_nan=False)
            stream.write('\n')
    print(json.dumps({'status': receipt['status'], 'receipt': str(receipt_path), 'error': receipt.get('error')}))
    return 0 if receipt['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
