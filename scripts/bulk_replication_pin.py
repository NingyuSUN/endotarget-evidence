#!/usr/bin/env python3
"""Fill the FILL_AT_FREEZE input pins of config/bulk_replication.json from the files on disk.

Run once before the freeze; refuses to change a pin that is already set to a different value.
``--root`` targets another project tree (tests). Standard library only.
"""
from pathlib import Path
import argparse, json, sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from followon_freeze import sha  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CONFIG = 'config/bulk_replication.json'


def main(root):
    root = Path(root); p = root / CONFIG; cfg = json.loads(p.read_text(encoding='utf-8'))
    changed = []
    for name, spec in cfg['inputs'].items():
        h = sha(root / spec['path'])
        if spec['sha256'] == 'FILL_AT_FREEZE': spec['sha256'] = h; changed.append(name)
        elif spec['sha256'] != h: sys.exit(f'pin already set and differs for {name}: {spec["sha256"][:12]} vs {h[:12]}')
    p.write_text(json.dumps(cfg, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'filled': changed}))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('--root'); a = ap.parse_args()
    main(Path(a.root).resolve() if a.root else ROOT)
