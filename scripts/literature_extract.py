#!/usr/bin/env python3
"""Convert downloaded JATS full-text XML into greppable text next to the XML, once.

For each reference/papers/PMID<id>/full.xml writes full.text.md (title, abstract, body by
section with headings, table captions and table cells as rows) and full.extract.json (XML
sha256, bytes, sections, characters, UTC). Reads only the local copies; no network. Standard
library only.
"""
from pathlib import Path
import hashlib, json, re, sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
PAPERS = ROOT / 'reference/papers'


def text_of(el):
    return re.sub(r'\s+', ' ', ''.join(el.itertext())).strip() if el is not None else ''


def strip_ns(root):
    for el in root.iter():
        if isinstance(el.tag, str) and '}' in el.tag: el.tag = el.tag.split('}', 1)[1]
    return root


def table_rows(tbl):
    rows = []
    for tr in tbl.iter('tr'):
        cells = [text_of(c) for c in tr if c.tag in ('td', 'th')]
        if any(cells): rows.append(' | '.join(cells))
    return rows


def convert(xml_path):
    data = xml_path.read_bytes(); root = strip_ns(ET.fromstring(data))
    art = root if root.tag == 'article' else root.find('.//article')
    out = []; sections = 0
    title = text_of(art.find('.//article-meta/title-group/article-title'))
    out.append(f'# {title}\n')
    ids = {i.get('pub-id-type'): text_of(i) for i in art.findall('.//article-meta/article-id')}
    out.append('IDs: ' + ', '.join(f'{k}={v}' for k, v in ids.items()) + '\n')
    for ab in art.findall('.//article-meta/abstract'):
        out.append('## Abstract\n'); out.append(text_of(ab) + '\n'); sections += 1
    body = art.find('body')
    def walk(sec, depth):
        nonlocal sections
        t = sec.find('title'); heading = text_of(t) if t is not None else ''
        if heading: out.append('#' * min(depth, 6) + ' ' + heading + '\n'); sections += 1
        for child in sec:
            if child.tag == 'sec': walk(child, depth + 1)
            elif child.tag == 'p': out.append(text_of(child) + '\n')
            elif child.tag in ('table-wrap', 'fig'):
                cap = text_of(child.find('caption')); lab = text_of(child.find('label'))
                out.append(f'[{child.tag}] {lab} {cap}\n')
                for tbl in child.iter('table'):
                    for r in table_rows(tbl): out.append('    ' + r)
                    out.append('')
            elif child.tag == 'list':
                for li in child.iter('list-item'): out.append('- ' + text_of(li))
    if body is not None:
        for child in body:
            if child.tag == 'sec': walk(child, 2)
            elif child.tag == 'p': out.append(text_of(child) + '\n')
    for tw in art.findall('.//floats-group/table-wrap'):
        out.append(f"[table-wrap floats] {text_of(tw.find('label'))} {text_of(tw.find('caption'))}\n")
        for tbl in tw.iter('table'):
            for r in table_rows(tbl): out.append('    ' + r)
    text = '\n'.join(out) + '\n'
    (xml_path.parent / 'full.text.md').write_text(text, encoding='utf-8')
    receipt = {'source': str(xml_path.relative_to(ROOT)).replace('\\', '/'), 'xml_sha256': hashlib.sha256(data).hexdigest(), 'xml_bytes': len(data),
               'title': title, 'sections': sections, 'characters': len(text), 'utc': datetime.now(timezone.utc).isoformat(), 'method': 'JATS XML, ElementTree'}
    (xml_path.parent / 'full.extract.json').write_text(json.dumps(receipt, indent=1), encoding='utf-8')
    return receipt


if __name__ == '__main__':
    targets = [Path(p).resolve() for p in sys.argv[1:]] or sorted(PAPERS.glob('PMID*/full.xml'))
    for x in targets:
        r = convert(x); print(f"{x.parent.name}: {r['characters']:>7} chars, {r['sections']:>3} sections | {r['title'][:70]}")
