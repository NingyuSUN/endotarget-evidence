#!/usr/bin/env python3
"""Generic freeze for a prespecified follow-on analysis described by a configuration file.

The configuration names the protocol, the files whose bytes are frozen (code, tests,
configuration), the inputs with their expected SHA-256 pins, and the output paths. The
freeze hashes the frozen files, verifies every input pin against the file on disk, records
the git head and refuses when any output or the results receipt already exists. The
analysis script later refuses to run if any frozen file changed. ``--root`` targets another
project tree (tests); ``--supersedes`` names an earlier receipt this one replaces.
"""
from pathlib import Path
import argparse, hashlib, json, subprocess, sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]


TEXT_SUFFIXES = ('.tsv', '.txt', '.csv', '.json', '.md', '.py', '.R', '.yml', '.yaml')


def sha(path):
    """SHA-256 of the file; text files are hashed with CRLF normalized to LF so that pins agree
    between a Windows checkout (autocrlf) and the Linux checkout that CI and the host use."""
    data = Path(path).read_bytes()
    if Path(path).suffix in TEXT_SUFFIXES: data = data.replace(b'\r\n', b'\n')
    return hashlib.sha256(data).hexdigest()


def main(root, config_rel, supersedes=None):
    root = Path(root); cfg = json.loads((root / config_rel).read_text(encoding='utf-8'))
    out = cfg['outputs']; phase = cfg['freeze']['phase']
    for rel in out['refuse_if_exists']:
        if (root / rel).exists(): sys.exit(f'refusing to freeze: output exists ({rel})')
    freeze_path = root / out['freeze_receipt']
    if freeze_path.exists() and supersedes != out['freeze_receipt']: sys.exit('freeze receipt exists; pass --supersedes to re-issue')
    files = list(dict.fromkeys(cfg['freeze']['files'] + [config_rel, cfg['protocol']]))
    frozen = {}
    for rel in files:
        if not (root / rel).exists():
            if cfg.get('synthetic'): continue
            sys.exit(f'frozen file missing: {rel}')
        frozen[rel] = sha(root / rel)
    pins = {}
    for name, spec in cfg['inputs'].items():
        p = root / spec['path']
        if not p.exists(): sys.exit(f'input missing: {spec["path"]}')
        if sha(p) != spec['sha256']: sys.exit(f'input hash differs from the configured pin: {spec["path"]}')
        pins[spec['path']] = spec['sha256']
    try: head = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=root, capture_output=True, text=True, check=True).stdout.strip()
    except Exception: head = None
    receipt = {'schema_version': 1, 'phase': phase, 'analysis': cfg['analysis'], 'created_utc': datetime.now(timezone.utc).isoformat(),
               'synthetic': bool(cfg.get('synthetic')), 'protocol': cfg['protocol'], 'configuration': config_rel, 'git_head': head,
               'frozen_sha256': frozen, 'pinned_inputs': pins, 'supersedes': supersedes,
               'attestations': cfg['freeze'].get('attestations', [])}
    freeze_path.parent.mkdir(parents=True, exist_ok=True)
    freeze_path.write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    return receipt


def verify(root, config_rel, analysis_script_rel):
    """Called by an analysis script before reading data: freeze present, phase right, frozen files unchanged, pins hold, no results yet."""
    root = Path(root); cfg = json.loads((root / config_rel).read_text(encoding='utf-8')); out = cfg['outputs']
    fp = root / out['freeze_receipt']
    if not fp.exists(): sys.exit('freeze receipt missing: run scripts/followon_freeze.py first')
    fr = json.loads(fp.read_text(encoding='utf-8'))
    if fr.get('phase') != cfg['freeze']['phase'] or fr.get('analysis') != cfg['analysis']: sys.exit('freeze receipt does not match this analysis')
    if bool(fr.get('synthetic')) != bool(cfg.get('synthetic')): sys.exit('synthetic flag differs between freeze and configuration')
    for rel, h in fr['frozen_sha256'].items():
        if sha(root / rel) != h: sys.exit(f'frozen file changed since the freeze: {rel}')
    if analysis_script_rel not in fr['frozen_sha256']: sys.exit('analysis script is not among the frozen files')
    for name, spec in cfg['inputs'].items():
        if sha(root / spec['path']) != spec['sha256'] or fr['pinned_inputs'].get(spec['path']) != spec['sha256']:
            sys.exit(f'input pin mismatch: {spec["path"]}')
    if (root / out['results_receipt']).exists(): sys.exit('results receipt exists; this analysis runs once')
    return cfg, fr, sha(fp)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', required=True, help='configuration path relative to the root')
    ap.add_argument('--root'); ap.add_argument('--supersedes')
    a = ap.parse_args()
    r = main(Path(a.root).resolve() if a.root else ROOT, a.config, a.supersedes)
    print(json.dumps({k: r[k] for k in ('phase', 'analysis', 'git_head')} | {'pinned': len(r['pinned_inputs']), 'frozen': len(r['frozen_sha256'])}))
