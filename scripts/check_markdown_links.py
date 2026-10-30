#!/usr/bin/env python3
"""Check that relative links in tracked Markdown files point at existing paths.

Standard library only. Historical copies under docs/archive/ keep the links of
their recorded version and are skipped. External URLs, mailto links and
in-page anchors are not checked. Exit status 1 when any link is broken.
``--root`` checks another git checkout (the boundary tests use a temporary one).
"""
import argparse, pathlib, re, subprocess, sys, urllib.parse

ROOT = pathlib.Path(__file__).resolve().parent.parent
SKIP_PREFIXES = ("docs/archive/",)
LINK = re.compile(r"\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")


def tracked_markdown(root):
    out = subprocess.run(["git", "ls-files", "*.md", "**/*.md"], cwd=root, check=True,
                         capture_output=True, text=True).stdout
    return sorted({line.strip() for line in out.splitlines() if line.strip()})


def main(root=ROOT):
    broken, checked = [], 0
    for rel in tracked_markdown(root):
        if rel.startswith(SKIP_PREFIXES):
            continue
        path = root / rel
        text = path.read_text(encoding="utf-8", errors="replace")
        for match in LINK.finditer(text):
            target = match.group(1)
            if target.startswith(("http://", "https://", "mailto:", "#")):
                continue
            target = urllib.parse.unquote(target.split("#", 1)[0])
            if not target:
                continue
            checked += 1
            if not (path.parent / target).exists():
                broken.append((rel, target))
    print(f"checked {checked} relative links; broken: {len(broken)}")
    for rel, target in broken:
        print(f"  {rel} -> {target}")
    return 1 if broken else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", help="git checkout to check (default: this repository)")
    args = ap.parse_args()
    sys.exit(main(pathlib.Path(args.root).resolve() if args.root else ROOT))
