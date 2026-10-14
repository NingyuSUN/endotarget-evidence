#!/usr/bin/env python3
"""Verify preserved science and the current local delivery, without running models.

Historical README hash exceptions apply only to two explicit old-path/old-hash
pairs. Current files are hashed directly. This is artifact and engineering
verification, not independent biological, causal, clinical or visual validation.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
from urllib.parse import unquote, urlsplit

from verify_local_migration import (
    HISTORICAL_ROOT, SOURCE_COMMIT, IntegrityChecker, git, require,
    verify_freezes, verify_delivery, verify_bootstrap, verify_lodo,
    verify_reported_results,
)
import check_local_container_nextflow as synthetic_checker
import check_local_primary_replay as primary_checker

OUTPUT = 'receipts/local_delivery_verification.json'
ARCHIVES = {
    'README.md': ('d73153d84255856c794640814c9c43a4b522c88b64e0aab8099dabea4588273e',
                  'docs/archive/README_srv01_4cc7580.md'),
    'README.zh-CN.md': ('f084154b11f8878422b8340ed70a017b431449263b9ed0474be23f80b48f9732',
                        'docs/archive/README.zh-CN_srv01_4cc7580.md'),
}
ARCHIVE_PATHS = frozenset(value[1] for value in ARCHIVES.values())
NEW_FIGURES = frozenset('reports/figures/' + name + '.svg' for name in (
    'donor_flow_v1', 'donor_composition_qc_v1', 'matched_null_sensitivity_v1', 'evidence_matrix_v1'))
CORE_CURRENT_FILES = NEW_FIGURES | ARCHIVE_PATHS | frozenset((
    'README.md', 'README.zh-CN.md', 'docs/METHOD_CARD.md', 'docs/LOCAL_EXECUTION.md',
    'docs/LOCAL_DELIVERY_ACCEPTANCE.md', 'docs/archive/INDEX.md',
    'reports/executive_brief_en.md', 'reports/executive_brief_en.pdf',
    'scripts/render_summary_figures_v1.py', 'scripts/render_executive_brief_pdf.py',
    'scripts/verify_local_delivery.py', 'tests/test_local_delivery.py',
    'receipts/summary_figures_v1.json', 'receipts/executive_brief_pdf.json',
    'receipts/local_container_build.json', 'receipts/local_container_manifest.json',
    'receipts/local_container_nextflow.json', 'receipts/local_primary_container_replay.json',
    'receipts/local_delivery_visual_review.json',
    'containers/statistics-local/Dockerfile', 'containers/statistics-local/.dockerignore',
    'containers/statistics-local/os-packages.lock.tsv',
))
CURRENT_SCOPES = frozenset(('scripts', 'tests', 'docs', 'reports', 'receipts',
                            'config', 'conf', 'containers', 'requirements', 'manifests'))


class DeliveryChecker(IntegrityChecker):
    def check_hash(self, name, expected, base=None):
        if base is None and isinstance(name, str) and name in ARCHIVES:
            old_hash, archive = ARCHIVES[name]
            require(expected == old_hash, 'Unexpected historical README hash: ' + name)
            return super().check_hash(archive, expected)
        return super().check_hash(name, expected, base)


def current_path(checker, name):
    """Exact POSIX spelling, selected delivery scopes, no traversal or symlinks."""
    require(isinstance(name, str) and bool(name) and '\x00' not in name,
            'Invalid current artifact path')
    p = PurePosixPath(name)
    require(not p.is_absolute() and p.as_posix() == name and '..' not in p.parts
            and '\\' not in name and ':' not in name,
            'Require normalized project-relative path: ' + name)
    require(name != OUTPUT, 'Delivery receipt cannot hash itself')
    require((len(p.parts) > 1 and p.parts[0] in CURRENT_SCOPES)
            or (len(p.parts) == 1 and not name.startswith('.') and p.suffix in ('.md', '.nf', '.config')),
            'Current artifact outside delivery scopes: ' + name)
    cursor = checker.root
    for part in p.parts:
        require(cursor.is_dir() and part in {entry.name for entry in cursor.iterdir()},
                'Missing artifact or incorrect path case: ' + name)
        cursor = cursor / part
        require(not cursor.is_symlink(), 'Current artifact symlink rejected: ' + name)
    require(cursor.resolve().is_relative_to(checker.root), 'Current artifact escapes project: ' + name)
    require(cursor.is_file(), 'Current artifact is not a file: ' + name)
    return cursor


def verify_current_files(checker, files):
    require(isinstance(files, list) and bool(files), 'Explicit current file list is required')
    require(len(files) == len(set(files)), 'Duplicate current file path')
    values = {name: checker.digest(current_path(checker, name)) for name in files}
    return {'count': len(values), 'sha256': values, 'self_receipt_excluded': OUTPUT}


def verify_required_current_files(checker, files):
    require(isinstance(files, list) and bool(files), 'Explicit current file list is required')
    missing = CORE_CURRENT_FILES - set(files)
    require(not missing, 'Missing core current delivery files: ' + ', '.join(sorted(missing)))
    result = verify_current_files(checker, files)
    result['required_core_files'] = sorted(CORE_CURRENT_FILES)
    return result


def verify_archives(checker):
    result = {}
    for name, (expected, archive) in ARCHIVES.items():
        original = subprocess.check_output(['git', '-C', str(checker.root), 'show', SOURCE_COMMIT + ':' + name])
        require(hashlib.sha256(original).hexdigest() == expected,
                'Source commit README identity differs: ' + name)
        current_path(checker, archive)
        checker.check_hash(name, expected)
        current = checker.digest(current_path(checker, name))
        result[name] = {'source_commit': SOURCE_COMMIT, 'archived_as': archive,
                        'historical_sha256': expected, 'current_sha256': current,
                        'root_readme_changed': current != expected}
    return {'readmes': result, 'policy': 'exact old relative path + old SHA256 + base=None only'}


def markdown_targets(text):
    # Code examples are not active links. Inline/reference Markdown links are.
    text = re.sub(r'^\s*(`{3,}|~{3,}).*?^\s*\1\s*$', '', text, flags=re.M | re.S)
    text = re.sub(r'`+[^`]*`+', '', text)
    targets = []
    for match in re.finditer(r'!?\[[^\]\n]*\]\(\s*(<[^>]+>|[^\s)]+)(?:\s+["\'][^\n]*?["\'])?\s*\)', text):
        targets.append(match.group(1).strip('<>'))
    definitions = {m.group(1).strip().lower(): m.group(2).strip('<>') for m in re.finditer(
        r'^\s{0,3}\[([^\]]+)\]:\s*(<[^>]+>|\S+)', text, flags=re.M)}
    for match in re.finditer(r'!?\[([^\]\n]+)\]\[([^\]\n]*)\]', text):
        key = (match.group(2) or match.group(1)).strip().lower()
        require(key in definitions, 'Undefined Markdown reference: ' + key)
        targets.append(definitions[key])
    for match in re.finditer(r'(?<!\])\[([^\]\n]+)\](?![\[(:])', text):
        key = match.group(1).strip().lower()
        if key in definitions:
            targets.append(definitions[key])
    return targets


def verify_markdown(checker, files):
    checked, skipped, local, planned = [], [], 0, []
    for name in files:
        if not name.endswith('.md'):
            continue
        p = current_path(checker, name)
        if name in ARCHIVE_PATHS:
            skipped.append(name)
            continue
        for target in markdown_targets(p.read_text(encoding='utf-8')):
            parts = urlsplit(target)
            if parts.scheme in ('http', 'https', 'mailto') or target.startswith('#'):
                continue
            require(not parts.scheme and not parts.netloc, 'Unsupported Markdown URL: ' + target)
            path = unquote(parts.path)
            require(bool(path) and '\\' not in path and '\x00' not in path,
                    'Invalid Markdown target: ' + target)
            selected = (p.parent / path).resolve()
            require(selected.is_relative_to(checker.root), 'Markdown target escapes project: ' + target)
            if selected == (checker.root / OUTPUT).resolve():
                planned.append({'document': name, 'target': OUTPUT})
                continue
            require(selected.exists(), 'Missing local Markdown target: ' + name + ' -> ' + target)
            local += 1
        checked.append(name)
    return {'files_checked': checked, 'local_targets_checked': local, 'explicit_archives_skipped': skipped,
            'planned_atomic_output_links': planned}


def verify_os_git_runtime(checker, build):
    """Compare the OS lock, build declaration and actual dpkg/Git inspection."""
    lock_name = 'containers/statistics-local/os-packages.lock.tsv'
    def packages(text, label, normalize_architecture=False):
        require(isinstance(text, str) and bool(text), 'Missing OS package table: ' + label)
        result = {}
        for row in text.splitlines():
            fields = row.split('\t')
            require(len(fields) == 2 and all(fields), 'Malformed OS package row: ' + label)
            name, version = fields
            if normalize_architecture and name.endswith(':amd64'):
                name = name[:-6]
            require(re.fullmatch(r'[a-z0-9][a-z0-9+.-]*', name) is not None
                    and version == version.strip() and not any(c.isspace() for c in version),
                    'Invalid OS package identity: ' + label)
            require(name not in result, 'Duplicate normalized OS package: ' + name)
            result[name] = version
        return result
    expected = packages(checker.resolve(lock_name).read_text(), 'lock')
    require(len(expected) == 37 and {'git', 'git-man', 'procps'}.issubset(expected),
            'Expected exactly 37 locked OS packages including Git/procps')
    declared = build.get('required_os_packages')
    require(isinstance(declared, dict) and declared == expected,
            'Declared OS package dictionary differs from 37-package lock')
    runtime = build.get('runtime_inspection', {})
    actual = packages(runtime.get('os_packages'), 'actual dpkg inspection', True)
    require(len(actual) == 37 and actual == expected,
            'Actual OS package versions differ from 37-package lock')
    require(build.get('provenance_git_available') is True,
            'Actual provenance Git availability has not been verified')
    version, head = runtime.get('git_version'), runtime.get('git_head')
    require(isinstance(version, str) and re.fullmatch(r'git version [0-9]+(?:\.[0-9]+)+[^\r\n]*', version.strip()) is not None,
            'Actual Git version is missing or invalid')
    locked_git_version = expected['git'].split(':')[-1].split('-', 1)[0]
    require(version.strip() == 'git version ' + locked_git_version,
            'Actual Git version differs from locked Git package')
    require(isinstance(head, str) and re.fullmatch(r'[0-9a-f]{40}', head) is not None,
            'Build inspection Git HEAD is missing or invalid')
    if 'source_git_head' in build:
        require(build['source_git_head'] == head, 'Host/container build Git HEAD binding differs')
    return {'installed_os_packages': 37, 'git_version': version.strip(),
            'git_head_at_build_inspection': head, 'provenance_git_available': True}


def verify_image_configuration(config, build, recipe):
    """Bind image defaults to the built recipe; task --user is checked separately."""
    require(config.get('architecture') == build.get('architecture') == 'amd64'
            and config.get('os') == build.get('os') == 'linux', 'Unexpected built platform')
    users = re.findall(r'^\s*USER\s+(\S+)\s*(?:#.*)?$', recipe, flags=re.I | re.M)
    default = config.get('config', {}).get('User')
    require(bool(users) and default == build.get('image_config', {}).get('User')
            == users[-1] == '57439:57439', 'Unexpected built image default user')
    return {'image_default_user': default, 'image_platform': 'linux/amd64',
            'task_runtime_user_checked_separately': '1003:1004'}


def verify_build_manifest(checker):
    build = checker.load('receipts/local_container_build.json')
    require(build.get('build_succeeded') is True and build.get('exit_code') == 0,
            'Local Docker build has not passed')
    require(build.get('scientific_lock_package_count') == 117
            and build.get('scientific_lock_exact_match') is True, 'Scientific package lock not verified')
    require('containers/statistics-local/os-packages.lock.tsv' in build['source_hashes'],
            'OS package lock is not bound to the built image')
    checker.check_map('local_build:source_hashes', build['source_hashes'])
    checker.check_hash(build['build_log'], build['build_log_sha256'])
    checker.check_hash('scripts/build_local_container.py', build['script_sha256'])
    lock = checker.resolve('requirements/statistical-runtime.explicit.txt').read_text().splitlines()
    expected = {line.split('#')[0].rsplit('/', 1)[1]: line.split('#')[1]
                for line in lock if line.startswith('https://')}
    installed = build['runtime_inspection']['scientific_packages']
    actual = {p['fn']: p['md5'] for p in installed}
    require(len(installed) == len(actual) == len(expected) == 117 and actual == expected,
            'Actual 117-package runtime inspection differs from lock')
    require(build.get('metrics_ps_available') is True
            and build['runtime_inspection']['ps_version'].startswith('ps from procps'),
            'Runtime metrics dependency not verified')
    os_git = verify_os_git_runtime(checker, build)
    head = os_git['git_head_at_build_inspection']
    require(git(checker.root, 'cat-file', '-t', head) == 'commit', 'Build inspection Git HEAD is not a commit')
    require(subprocess.run(['git', '-C', str(checker.root), 'merge-base', '--is-ancestor', head, 'HEAD'],
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE).returncode == 0,
            'Build inspection Git HEAD is not an ancestor of current HEAD')
    manifest = checker.load('receipts/local_container_manifest.json')
    require(manifest.get('verified') is True and manifest.get('registry_stopped') is True,
            'Local manifest verification/registry shutdown has not passed')
    require(manifest.get('external_publication') is False, 'Unexpected external publication scope')
    checker.check_map('local_manifest:artifacts', manifest['artifacts'])
    checker.check_hash('scripts/verify_local_container_manifest.py', manifest['script_sha256'])
    for field in ('push_manifest', 'platform_manifest'):
        value = manifest[field]
        digest = value['requested_digest']
        require(re.fullmatch(r'sha256:[0-9a-f]{64}', digest) is not None
                and value['header_digest'] == value['raw_bytes_sha256'] == digest,
                'Manifest digest/header identity mismatch')
        checker.check_hash(value['file'], digest.split(':')[1])
    index = checker.load(manifest['push_manifest']['file'])
    platform = checker.load(manifest['platform_manifest']['file'])
    candidates = [m for m in index.get('manifests', []) if m.get('platform', {}).get('os') == 'linux'
                  and m.get('platform', {}).get('architecture') == 'amd64']
    require(len(candidates) == 1 and candidates[0]['digest'] == manifest['platform_digest']
            == manifest['platform_manifest']['requested_digest'], 'Selected Linux amd64 manifest differs')
    require(platform['config']['digest'] == manifest['image_config_digest'] == manifest['manifest_config_digest'],
            'Manifest config chain differs')
    config_files = [name for name in manifest['artifacts'] if Path(name).name == 'image-config.json']
    require(len(config_files) == 1, 'Expected one image-config blob')
    checker.check_hash(config_files[0], manifest['image_config_digest'].split(':')[1])
    config = checker.load(config_files[0])
    recipe_name = 'containers/statistics-local/Dockerfile'
    require(recipe_name in build['source_hashes'], 'Built recipe is not bound to image source hashes')
    image_configuration = verify_image_configuration(config, build, checker.resolve(recipe_name).read_text())
    require(build['image_id'] == manifest['image_id'] == manifest['build_image_id']
            and build['image_id'] in (manifest['push_digest'], manifest['platform_digest'], manifest['image_config_digest']),
            'Build/manifest image identity differs')
    image = manifest['image_ref']
    synthetic_checker.validate_image(image)
    require(image.endswith('@' + manifest['push_digest']) and image in manifest['repo_digests'],
            'Requested image is not the verified registry manifest')
    return {'image': image, 'image_id': build['image_id'], 'installed_scientific_packages': 117,
            **os_git, **image_configuration, 'registry_stopped': True, 'build_receipt_sha256': checker.digest('receipts/local_container_build.json')}


def verify_protected_inputs(checker, receipt, label):
    before, after = receipt['input_sha256_before'], receipt['input_sha256_after']
    require(before == after and bool(before), label + ' protected inputs changed')
    require(receipt.get('original_inputs_and_native_outputs_preserved') is True,
            label + ' did not preserve original evidence')
    checker.check_map(label + ':protected_inputs', before)


def verify_containers(checker, receipt, taskdirs):
    require(len(taskdirs) == len(set(taskdirs)) == 5, 'Expected exactly five distinct fresh tasks')
    items = receipt['container_inspections']
    require(len(items) == 5 and not receipt.get('container_observer_errors'), 'Actual container coverage failed')
    found = set()
    for item in items:
        value = item['inspect']
        require(value['Config']['Image'] == receipt['image'] and value['Image'] == receipt['image_id'],
                'Actual container image identity differs')
        matches = set(taskdirs).intersection(m['Destination'] for m in value.get('Mounts', []))
        require(len(matches) == 1, 'Container must match exactly one task workdir')
        work = next(iter(matches))
        checker.resolve(work)
        synthetic_checker.check_container_inspect(value, HISTORICAL_ROOT, Path(work))
        require(work not in found, 'Duplicate actual task container')
        found.add(work)
    require(found == set(taskdirs), 'Missing actual task container')


def verify_nextflow(checker, build_manifest):
    obj = checker.load('receipts/local_container_nextflow.json')
    require(obj.get('status') == 'PASS' and obj.get('docker_execution_verified') is True
            and obj.get('fresh_tasks_completed') == 5 and obj.get('resumed_tasks_cached') == 5,
            'Actual synthetic Nextflow fresh/cache gate has not passed')
    require(obj.get('discovery_tsvs_compared') == 24 and obj.get('replication_tsvs_compared') == 10,
            'Synthetic comparison coverage differs')
    require(obj['image'] == build_manifest['image'] and obj['image_id'] == build_manifest['image_id'],
            'Synthetic run does not use verified image')
    checker.check_hash('receipts/local_container_build.json', obj['build_receipt_sha256'])
    verify_protected_inputs(checker, obj, 'synthetic')
    expected = {('discovery', False), ('discovery', True), ('replication', False), ('replication', True)}
    require(len(obj['launches']) == 4 and {(r['pipeline'], r['resume']) for r in obj['launches']} == expected,
            'Synthetic fresh/cache pipeline coverage differs')
    taskdirs = []
    for run in obj['launches']:
        require(run['exit_code'] == 0, 'Synthetic launch failed')
        _, _, stages, count = synthetic_checker.FIXTURES[run['pipeline']]
        rows = synthetic_checker.check_trace(checker.resolve(run['trace']), stages, run['resume'])
        require(rows == run['trace_rows'], 'Saved synthetic trace rows differ')
        fixture = synthetic_checker.FIXTURES[run['pipeline']][0]
        comparison = synthetic_checker.compare_tree(checker.resolve('work/synthetic/' + fixture + '/results_first'),
                                                     checker.resolve(run['output']), count)
        require(comparison == run['all_scientific_tsv_comparisons'], 'Synthetic compared table provenance differs')
        for task in run['tasks']:
            work = task['workdir']
            checker.check_hash(work + '/.command.run', task['command_run_sha256'])
            checker.check_hash(work + '/.command.sh', task['command_sh_sha256'])
            if not run['resume']:
                taskdirs.append(work)
        if run['resume']:
            require(run['cached_output_bytes_identical'] is True, 'Cached synthetic outputs differ')
    verify_containers(checker, obj, taskdirs)
    return {'fresh_tasks': 5, 'cached_tasks': 5, 'scientific_tsvs_per_fresh_or_cached_comparison': 34}


def verify_primary(checker, build_manifest):
    obj = checker.load('receipts/local_primary_container_replay.json')
    require(obj.get('status') == 'PASS' and obj.get('docker_execution_verified') is True
            and obj.get('real_primary_outputs_reproduced') is True and obj.get('fresh_tasks_completed') == 5,
            'Actual frozen primary container replay has not passed')
    require(obj['image'] == build_manifest['image'] and obj['image_id'] == build_manifest['image_id'],
            'Primary run does not use verified image')
    require(obj.get('discovery_tsvs_compared') == 24 and obj.get('HECA_tsvs_compared') == 10,
            'Primary comparison coverage differs')
    for field, path in [('build_receipt_sha256', 'receipts/local_container_build.json'),
                        ('synthetic_receipt_sha256', 'receipts/local_container_nextflow.json')]:
        checker.check_hash(path, obj[field])
    require(obj.get('new_scientific_validation') is False and obj.get('bootstrap_500_rerun') is False
            and obj.get('LODO_22_rerun') is False and obj.get('raw_import_rerun') is False,
            'Primary replay scope changed')
    verify_protected_inputs(checker, obj, 'primary')
    require(len(obj['launches']) == 2 and {r['pipeline'] for r in obj['launches']} == set(primary_checker.JOBS),
            'Primary pipeline coverage differs')
    taskdirs = []
    for run in obj['launches']:
        job = primary_checker.JOBS[run['pipeline']]
        require(run['source'] == job['source'] and run['resume'] is False and run['exit_code'] == 0,
                'Primary launch/source scope differs')
        trace = obj['run_directory'] + '/' + run['pipeline'] + '_trace.tsv'
        require(synthetic_checker.check_trace(checker.resolve(trace), job['stages'], False) == run['trace_rows'],
                'Saved primary trace rows differ')
        comparison = primary_checker.compare_outputs(checker.resolve(job['source']), checker.resolve(run['output']), job)
        require(all(comparison[key] == run[key] for key in comparison), 'Primary output comparison provenance differs')
        for task in run['tasks']:
            work = task['workdir']
            checker.check_hash(work + '/.command.run', task['command_run_sha256'])
            checker.check_hash(work + '/.command.sh', task['command_sh_sha256'])
            taskdirs.append(work)
    verify_containers(checker, obj, taskdirs)
    return {'fresh_tasks': 5, 'scientific_tsvs_compared': 34, 'real_primary_outputs_reproduced': True,
            'new_scientific_validation': False, 'bootstrap_or_LODO_reexecuted': False}


def verify_figures(checker):
    obj = checker.load('receipts/summary_figures_v1.json')
    require(set(obj['output_sha256']) == NEW_FIGURES and obj.get('scientific_fits_executed') == 0
            and obj.get('new_disease_tests') == 0, 'Figure coverage/render-only scope differs')
    checker.check_map('summary_figures:inputs', obj['input_sha256'])
    checker.check_map('summary_figures:outputs', obj['output_sha256'])
    checker.check_hash('scripts/render_summary_figures_v1.py', obj['renderer_sha256'])
    return {'new_figures': sorted(NEW_FIGURES), 'scientific_fits_executed': 0}


def verify_pdf(checker):
    obj = checker.load('receipts/executive_brief_pdf.json')
    source, output = obj['source'].replace('\\', '/'), obj['output'].replace('\\', '/')
    require(source == 'reports/executive_brief_en.md' and output == 'reports/executive_brief_en.pdf',
            'Unexpected English brief source/output')
    require(type(obj.get('pages')) is int and obj['pages'] == 2 and obj.get('new_scientific_analysis') is False,
            'PDF receipt is not a two-page render-only delivery')
    checker.check_hash(source, obj['source_sha256'])
    checker.check_hash(output, obj['output_sha256'])
    checker.check_hash('scripts/render_executive_brief_pdf.py', obj['script_sha256'])
    pdf = checker.resolve(output).read_bytes()
    require(pdf.startswith(b'%PDF-') and len(re.findall(rb'/Type\s*/Page\b', pdf)) == 2,
            'Expected exactly two explicit ReportLab PDF page objects')
    return {'pages': 2, 'source': source, 'output': output, 'visual_review_checked_separately': True}


def verify_visual_review(checker, filename):
    required = set(NEW_FIGURES) | {'reports/executive_brief_en.pdf'}
    current_path(checker, filename)
    obj = checker.load(filename)
    require(obj.get('status') == 'PASS' and obj.get('reviewer') == 'Codex visual inspection'
            and obj.get('pdf_pages') == 2, 'Visual inspection receipt has not passed')
    require(set(obj.get('artifacts', {})) == required, 'Visual inspection artifact coverage differs')
    checker.check_map('visual_review:artifacts', obj['artifacts'])
    if obj.get('previews'):
        checker.check_map('visual_review:previews', obj['previews'])
    return {'receipt': filename, 'artifacts': sorted(required), 'automated_visual_proof': False,
            'scope': 'Integrity of recorded human/model visual inspection; no automated proof of appearance'}


def run(root, files, visual_audit):
    checker = DeliveryChecker(root)
    report = {'schema_version': 1, 'created_utc': datetime.now(timezone.utc).isoformat(),
              'source_commit': SOURCE_COMMIT, 'local_root': str(checker.root),
              'local_head_before_receipt': git(checker.root, 'rev-parse', 'HEAD'),
              'scientific_validation': False, 'therapeutic_validation': False,
              'new_scientific_fits_executed': 0, 'checks': {}, 'errors': []}
    subprocess.check_call(['git', '-C', str(checker.root), 'merge-base', '--is-ancestor',
                           SOURCE_COMMIT, report['local_head_before_receipt']])
    stages = [('archives', verify_archives), ('freezes', verify_freezes), ('historical_delivery', verify_delivery),
              ('bootstrap', verify_bootstrap), ('LODO', verify_lodo), ('reported_results', verify_reported_results),
              ('build_manifest', verify_build_manifest), ('figures', verify_figures), ('two_page_brief', verify_pdf),
              ('current_files', lambda c: verify_required_current_files(c, files)),
              ('active_markdown', lambda c: verify_markdown(c, files)),
              ('visual_review', lambda c: verify_visual_review(c, visual_audit))]
    for name, function in stages:
        try:
            report['checks'][name] = {'status': 'PASS', **function(checker)}
        except (ValueError, OSError, KeyError, TypeError, IndexError, subprocess.CalledProcessError) as error:
            report['checks'][name] = {'status': 'FAIL'}
            report['errors'].append({'stage': name, 'error': str(error)})
    binding = report['checks'].get('build_manifest', {})
    for name, function in [('synthetic_container', verify_nextflow), ('primary_container', verify_primary)]:
        try:
            require(binding.get('status') == 'PASS', 'Build/manifest prerequisite failed')
            report['checks'][name] = {'status': 'PASS', **function(checker, binding)}
        except (ValueError, OSError, KeyError, TypeError, IndexError, subprocess.CalledProcessError) as error:
            report['checks'][name] = {'status': 'FAIL'}
            report['errors'].append({'stage': name, 'error': str(error)})
    verifier_sha256 = checker.digest(Path(__file__).resolve())
    report.update(status='PASS_LOCAL_DELIVERY_INTEGRITY' if not report['errors'] else 'FAIL_LOCAL_DELIVERY_INTEGRITY',
                  hash_checks=checker.hash_checks, unique_files_hashed=checker.unique_files_hashed,
                  unique_bytes_hashed=sum(x[2] for x in checker.cache.values()), completed_hash_maps=checker.map_checks,
                  verifier_sha256=verifier_sha256)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--files', nargs='+', required=True, help='Exact project-relative current file list; no globbing')
    parser.add_argument('--visual-audit', default='receipts/local_delivery_visual_review.json')
    args = parser.parse_args()
    report = run(args.root, args.files, args.visual_audit)
    output = args.root.resolve() / OUTPUT
    temporary = output.with_suffix('.json.partial')
    temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    temporary.replace(output)
    print(json.dumps({key: report[key] for key in ('status', 'hash_checks', 'unique_files_hashed', 'errors')}))
    return 0 if not report['errors'] else 1


if __name__ == '__main__':
    sys.exit(main())
