#!/usr/bin/env python3
"""Freeze the two-cohort combined bound analysis before any value is computed.

Hashes the protocol, configuration, analysis and freeze scripts and tests; verifies that
the configured input pins match the files and the MDE results receipt; records the git
head; refuses when any combined output or results receipt already exists. Writes
receipts/combined_bound_freeze.json. ``--root`` targets another project tree (tests);
``--supersedes`` names an earlier freeze receipt that this one replaces.
"""
from pathlib import Path
import argparse, hashlib, json, subprocess, sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
PHASE = 'COMBINED_BOUND_RULES_FROZEN_PRE_ANALYSIS'
CODE = ['config/combined_bound.json', 'scripts/combined_bound_analysis.py', 'scripts/combined_bound_freeze.py', 'tests/test_combined_bound.py']


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(root, supersedes=None):
    root = Path(root); cfg = json.loads((root / 'config/combined_bound.json').read_text(encoding='utf-8'))
    out = cfg['outputs']
    for rel in (out['results_receipt'], out['per_gene'], out['summary']):
        if (root / rel).exists(): sys.exit(f'refusing to freeze: output exists ({rel})')
    freeze_path = root / out['freeze_receipt']
    if freeze_path.exists() and supersedes != out['freeze_receipt']: sys.exit('freeze receipt exists; pass --supersedes to re-issue')
    frozen = {rel: sha(root / rel) for rel in CODE + [cfg['protocol']] if (root / rel).exists()}
    missing = [rel for rel in CODE + [cfg['protocol']] if rel not in frozen]
    if missing and not cfg.get('synthetic'): sys.exit(f'frozen files missing: {missing}')
    pins = {}
    for c in ('discovery', 'heca'):
        spec = cfg['inputs'][c]; p = root / spec['path']
        if not p.exists(): sys.exit(f'input missing: {spec["path"]}')
        if sha(p) != spec['sha256']: sys.exit(f'input hash differs from the configured pin: {spec["path"]}')
        pins[spec['path']] = spec['sha256']
    receipt_spec = cfg['inputs'].get('pinning_receipt')
    pin_check = None
    if receipt_spec and receipt_spec.get('must_list_inputs') and not cfg.get('synthetic'):
        rec = json.loads((root / receipt_spec['path']).read_text(encoding='utf-8'))
        listed = rec.get('artifacts', {})
        for rel, h in pins.items():
            if listed.get(rel) != h: sys.exit(f'MDE results receipt does not list {rel} with the configured hash')
        pin_check = {'receipt': receipt_spec['path'], 'receipt_sha256': sha(root / receipt_spec['path']), 'inputs_listed_with_matching_hash': True}
    try:
        head = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=root, capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        head = None
    receipt = {'schema_version': 1, 'phase': PHASE, 'created_utc': datetime.now(timezone.utc).isoformat(), 'synthetic': bool(cfg.get('synthetic')),
               'protocol': cfg['protocol'], 'git_head': head, 'frozen_sha256': frozen, 'pinned_inputs': pins, 'pinning_receipt_check': pin_check,
               'supersedes': supersedes,
               'attestations': ['no combined estimate, p value, count or minimum detectable effect had been computed when this receipt was written',
                                'the analysis reads only the two pinned per-gene tables; no count matrix, design or model is touched',
                                'labels and sensitivities are those of the protocol; sensitivities never change a label']}
    freeze_path.parent.mkdir(parents=True, exist_ok=True)
    freeze_path.write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    return receipt


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root'); ap.add_argument('--supersedes')
    a = ap.parse_args()
    r = main(Path(a.root).resolve() if a.root else ROOT, a.supersedes)
    print(json.dumps({k: r[k] for k in ('phase', 'git_head', 'pinned_inputs')}))
