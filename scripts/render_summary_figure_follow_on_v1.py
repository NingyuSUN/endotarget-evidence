#!/usr/bin/env python3
"""Render one summary figure of the follow-on analyses from published tables and receipts.

Standard library only. Reads the published minimum-detectable-effect summaries and
donor-subsampling curves, the secretory-program summaries, the positive-control
evaluation receipt and the v1 dossier pool; cross-checks them; writes an SVG and a
receipt with input, script and output hashes plus every number drawn. No model is
fitted or rerun, and no scientific value is changed. ``--root`` points the script
at another project tree (the boundary tests use a temporary copy).
"""
from pathlib import Path
import argparse, csv, hashlib, html, json
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
COHORTS = ('discovery', 'heca')
INPUTS = [f'reports/mde_analysis/v1_{c}_summary.json' for c in COHORTS] \
    + [f'reports/mde_analysis/phase_{c}_summary.json' for c in COHORTS] \
    + [f'reports/mde_analysis/curve_{c}_curve_summary.tsv' for c in COHORTS] \
    + [f'reports/secretory_program/{c}_summary.json' for c in COHORTS] \
    + ['receipts/phase_positive_control_evaluation.json', 'reports/replication_v2_dossier_pool.tsv']
OUTPUT = 'reports/figures/follow_on_summary_v1.svg'
RECEIPT = 'receipts/summary_figure_follow_on_v1.json'
NAME = {'discovery': 'Discovery', 'heca': 'HECA'}
COLOR = {'discovery': '#4477AA', 'heca': '#EE7733'}
TINT = {'discovery': '#D5DFEA', 'heca': '#FBE3D6'}   # solid tints instead of opacity, for every renderer
GREY = '#9AA5B1'; INK = '#243247'; MUTED = '#5B6B7C'; GRID = '#dce2e9'


def jload(root, rel): return json.loads((root / rel).read_text(encoding='utf-8'))
def tsv(root, rel):
    with (root / rel).open(newline='', encoding='utf-8') as f: return list(csv.DictReader(f, delimiter='\t'))
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def esc(s): return html.escape(str(s), quote=True)
def signed(x, digits=2): return f'{x:+.{digits}f}'.replace('-', '−')
def text(x, y, s, size=13, fill=INK, extra=''):
    return f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" {extra}>{esc(s)}</text>'
def rect(x, y, w, h, fill, extra=''):
    return f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" fill="{fill}" {extra}/>'
def line(x1, y1, x2, y2, color=GRID, width=1, extra=''):
    return f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{color}" stroke-width="{width}" {extra}/>'
def circle(x, y, r, fill): return f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r}" fill="{fill}"/>'
def poly(points, fill='none', stroke='none', width=1, closed=False):
    tag = 'polygon' if closed else 'polyline'
    pts = ' '.join(f'{x:.2f},{y:.2f}' for x, y in points)
    return f'<{tag} points="{pts}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>'


def load_and_check(root):
    """Read every input and cross-check the numbers before anything is drawn."""
    v1 = {c: jload(root, f'reports/mde_analysis/v1_{c}_summary.json') for c in COHORTS}
    phase = {c: jload(root, f'reports/mde_analysis/phase_{c}_summary.json') for c in COHORTS}
    curve = {c: tsv(root, f'reports/mde_analysis/curve_{c}_curve_summary.tsv') for c in COHORTS}
    sp = {c: jload(root, f'reports/secretory_program/{c}_summary.json') for c in COHORTS}
    ev = jload(root, 'receipts/phase_positive_control_evaluation.json')
    pool = tsv(root, 'reports/replication_v2_dossier_pool.tsv')
    for c in COHORTS:
        assert v1[c]['synthetic'] is False and phase[c]['synthetic'] is False and sp[c]['synthetic'] is False, 'synthetic input'
        assert v1[c]['coefficient'] == 'ENDO_minus_non_ENDO' and sp[c]['M1_primary']['coefficient'] == 'SEC_x_ENDO', 'wrong coefficient'
    disease_sig = {c: v1[c]['observed']['genes_BH_q_lt_0_05'] for c in COHORTS}
    phase_sig = {c: phase[c]['observed']['genes_BH_q_lt_0_05'] for c in COHORTS}
    P1, P2 = ev['endpoints']['P1'], ev['endpoints']['P2']
    assert disease_sig == {'discovery': 0, 'heca': 0}, 'v1 disease fits are expected to have no gene at q < 0.05'
    assert phase_sig['discovery'] == P1['discovery_BH_q_lt_0_05'] and phase_sig['heca'] == P1['heca_BH_q_lt_0_05'], \
        'phase gene counts differ between the MDE summaries and the evaluation receipt'
    assert P2['selected'] == 20 and P2['pass'] is True and P1['pass'] is True, 'positive-control endpoints not as recorded'
    b0 = [r for r in pool if r['in_B0_top20'] == 'TRUE']
    assert len(b0) == 20 and len({r['gene_id'] for r in b0}) == 20, 'frozen B0 top-20 must be exactly 20 distinct genes'
    v1_transfer = {'same_sign': sum(r['agree'] == 'TRUE' for r in b0),
                   'evaluable': sum(r['evaluable'] == 'TRUE' for r in b0),
                   'same_sign_and_external_q_lt_0_05': sum(r['agree'] == 'TRUE' and r['external_q'] not in ('', 'NA')
                                                           and float(r['external_q']) < 0.05 for r in b0)}
    assert v1_transfer['evaluable'] == 20 and v1_transfer['same_sign_and_external_q_lt_0_05'] == 0, 'v1 transfer counts not as recorded'
    mde = {c: {s: v1[c]['scenarios'][s]['quantiles_fold']['q0.5'] for s in ('alpha_0.05', 'bh_k100', 'bh_k1')} for c in COHORTS}
    family = {c: v1[c]['family_size'] for c in COHORTS}
    curves = {}
    for c in COHORTS:
        rows = [r for r in curve[c] if r['design'] == 'primary']
        assert rows and all(r['cohort'] == c for r in rows), 'curve rows missing or mislabeled'
        curves[c] = [{'n': int(r['level_n_per_arm']), 'valid': int(r['valid']), 'draws': int(r['draws']),
                      'median': None if r['recall_full_median'] == 'NA' else float(r['recall_full_median']),
                      'q025': None if r['recall_full_q025'] == 'NA' else float(r['recall_full_q025']),
                      'q975': None if r['recall_full_q975'] == 'NA' else float(r['recall_full_q975'])} for r in rows]
        assert all((p['valid'] == 0) == (p['median'] is None) for p in curves[c]), 'recall present without valid draws'
    theta = {}
    for c in COHORTS:
        m, g = sp[c]['M1_primary'], sp[c]['validity_gate']
        assert m['status'] == 'ESTIMATED' and g['pass'] is True and sp[c]['label'] == 'no_detectable_difference', \
            'secretory program result not as recorded'
        d = sp[c]['donors']['by_phase_disease']
        theta[c] = {'estimate': m['estimate'], 'ci_low': m['ci_low'], 'ci_high': m['ci_high'], 'p': m['p'],
                    'mde': m['mde_log'], 'auc': g['auc'], 'shift': g['model']['estimate'],
                    'secretory_endo': d['secretory ENDO'], 'secretory_non_endo': d['secretory non-ENDO'],
                    'n_donors': sp[c]['donors']['n'],
                    'max_attenuation_share': -m['ci_low'] / g['model']['estimate'],
                    'max_enhancement_share': m['ci_high'] / g['model']['estimate']}
    at11 = {c: next(p for p in curves[c] if p['n'] == 11) for c in COHORTS}
    return {'disease_sig': disease_sig, 'phase_sig': phase_sig, 'P2': P2, 'v1_transfer': v1_transfer, 'mde': mde,
            'family': family, 'curves': curves, 'theta': theta, 'at11': at11}


def draw(d):
    disease_sig, phase_sig, P2, v1_transfer = d['disease_sig'], d['phase_sig'], d['P2'], d['v1_transfer']
    mde, curves, theta, at11 = d['mde'], d['curves'], d['theta'], d['at11']
    W, H = 1500, 650
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-labelledby="title desc">',
             f'<title id="title">{esc("What the null means: follow-on analyses of the v1 benchmark")}</title>',
             f'<desc id="desc">{esc("Three panels rendered from published tables: detection and transfer of the disease and cycle-phase contrasts, recovery of the phase program by donor number, and the secretory program disease interaction with its detectable effect.")}</desc>',
             '<rect width="100%" height="100%" fill="white"/><g font-family="DejaVu Sans, Arial, sans-serif">']
    parts += [text(35, 40, 'What the null means: the v1 pipeline under four prespecified follow-on analyses', 24),
              text(35, 66, 'Rendered from published summaries and receipts; no model rerun. Disease contrast: 22 proliferative donors per cohort. Follow-on contrasts add the secretory arms: 35 (discovery) and 42 (HECA) donors.', 14, MUTED)]
    LINES = (486, 506, 526, 556)   # text block under each panel

    # Panel A: detection and transfer
    ax, top, bottom = 35, 160, 430
    left, right = 100, 470
    ymax = 8000
    parts += [text(ax, 118, 'A. Genes at BH q < 0.05: disease versus cycle phase', 17),
              text(ax, 138, 'Same estimator, filter, normalization and covariates; only the contrast differs', 12, MUTED)]
    for v in range(0, ymax + 1, 2000):
        y = bottom - v / ymax * (bottom - top)
        parts += [line(left, y, right, y), text(left - 10, y + 4, f'{v:,}', 11, MUTED, 'text-anchor="end"')]
    for i, c in enumerate(COHORTS):
        cx = left + 95 + i * 180
        for j, (val, col, lab) in enumerate([(disease_sig[c], GREY, 'disease'), (phase_sig[c], COLOR[c], 'phase')]):
            x = cx - 60 + j * 65; h = val / ymax * (bottom - top)
            if val == 0:
                parts += [line(x, bottom - 1.5, x + 55, bottom - 1.5, col, 3)]
            else:
                parts += [rect(x, bottom - h, 55, h, col)]
            parts += [text(x + 27.5, bottom - h - 8, f'{val:,}', 12, INK, 'text-anchor="middle"'),
                      text(x + 27.5, bottom + 16, lab, 11, MUTED, 'text-anchor="middle"')]
        unit = 'scRNA' if c == 'discovery' else 'snRNA'
        acc = 'GSE266265' if c == 'discovery' else 'E-MTAB-14039'
        parts += [text(cx + 2, bottom + 36, f'{NAME[c]} ({acc}, {unit})', 12, INK, 'text-anchor="middle"')]
    parts += [text(ax, LINES[0], 'Frozen expression (B0) top-20 transferred to HECA:', 12.5),
              text(ax, LINES[1], f"disease: {v1_transfer['same_sign']} / {v1_transfer['evaluable']} / 20 same direction / evaluable / selected; {v1_transfer['same_sign_and_external_q_lt_0_05']} at external q < 0.05", 11.5, MUTED),
              text(ax, LINES[2], f"phase: {P2['same_sign']} / {P2['evaluable']} / {P2['selected']} same direction; {P2['same_sign_and_external_q_lt_0_05']} at external q < 0.05", 11.5, MUTED),
              text(ax, LINES[3], 'Reading: the null is not pipeline insensitivity; phase genes are not a shortlist.', 11.5, INK)]

    # Panel B: subsampling recall curve
    bx = 520; left, right = 585, 975; top, bottom = 160, 430
    xmin, xmax = 4, 21
    def X(n): return left + (n - xmin) / (xmax - xmin) * (right - left)
    def Y(r): return bottom - r * (bottom - top)
    parts += [text(bx, 118, 'B. Recovery of the phase program by donor number', 17),
              text(bx, 138, 'Subsampled studies of n + n donors re-analyzed with the v1 estimator (primary design)', 12, MUTED)]
    for r in (0, 0.25, 0.5, 0.75, 1.0):
        parts += [line(left, Y(r), right, Y(r)), text(left - 10, Y(r) + 4, f'{int(r * 100)}%', 11, MUTED, 'text-anchor="end"')]
    for n in (5, 6, 8, 10, 11, 13, 16, 20):
        parts += [text(X(n), bottom + 16, str(n), 11, MUTED, 'text-anchor="middle"')]
    parts += [text((left + right) / 2, bottom + 36, 'donors per arm (secretory and proliferative)', 12, INK, 'text-anchor="middle"'),
              text(left - 48, (top + bottom) / 2, 'recall', 12, INK, f'transform="rotate(-90 {left - 48:.1f} {(top + bottom) / 2:.1f})" text-anchor="middle"')]
    for c in COHORTS:   # tinted bands first, so lines and points stay visible
        pts = [p for p in curves[c] if p['median'] is not None]
        band = [(X(p['n']), Y(p['q975'])) for p in pts] + [(X(p['n']), Y(p['q025'])) for p in reversed(pts)]
        parts += [poly(band, fill=TINT[c], closed=True)]
    parts += [line(X(11), top, X(11), bottom, '#5B6B7C', 1, 'stroke-dasharray="4 4"'),
              text(X(11) + 5, top + 12, '11 + 11 = a 22-donor study', 11, MUTED)]
    for c in COHORTS:
        pts = [p for p in curves[c] if p['median'] is not None]
        parts += [poly([(X(p['n']), Y(p['median'])) for p in pts], stroke=COLOR[c], width=2.2)]
        parts += [circle(X(p['n']), Y(p['median']), 3.5, COLOR[c]) for p in pts]
        for p in curves[c]:
            if p['median'] is None:
                parts += [text(X(p['n']), bottom - 6 - (12 if c == 'heca' else 0), '×', 13, COLOR[c], 'text-anchor="middle"')]
        m = at11[c]
        parts += [text(X(11) + 7, Y(m['median']) - 10, f"{NAME[c]} {m['median'] * 100:.0f}% [{m['q025'] * 100:.0f}–{m['q975'] * 100:.0f}%]", 11.5, COLOR[c])]
    lx, ly = X(16.2), bottom - 52
    parts += [rect(lx, ly, 12, 12, COLOR['discovery']), text(lx + 17, ly + 11, 'Discovery (≤13 per arm)', 11, MUTED),
              rect(lx, ly + 20, 12, 12, COLOR['heca']), text(lx + 17, ly + 31, 'HECA (≤20 per arm)', 11, MUTED),
              text(bx, LINES[0], 'Median and 2.5th–97.5th percentiles over 100 draws per level; × = design not estimable.', 11.5, MUTED),
              text(bx, LINES[1], 'Per-gene detectable effect of the v1 disease fits (median gene, discovery / HECA):', 12.5),
              text(bx, LINES[2], f"{mde['discovery']['alpha_0.05']:.1f}× / {mde['heca']['alpha_0.05']:.1f}× (α 0.05) · {mde['discovery']['bh_k100']:.1f}× / {mde['heca']['bh_k100']:.1f}× (100 BH discoveries) · {mde['discovery']['bh_k1']:.1f}× / {mde['heca']['bh_k1']:.1f}× (one)", 11.5, MUTED),
              text(bx, LINES[3], 'Reading: 22 donors recover even a large known program inconsistently.', 11.5, INK)]

    # Panel C: secretory program interaction
    cx0 = 1025; left, right = 1085, 1450; top, bottom = 185, 360
    tmin, tmax = -0.75, 0.45
    def TX(t): return left + (t - tmin) / (tmax - tmin) * (right - left)
    parts += [text(cx0, 118, 'C. Secretory stromal program: disease interaction θ', 17),
              text(cx0, 138, 'Scored on the other cohort’s frozen phase program; 1 = the full secretory shift', 12, MUTED)]
    for t in (-0.6, -0.4, -0.2, 0.0, 0.2, 0.4):
        parts += [line(TX(t), top, TX(t), bottom), text(TX(t), bottom + 16, signed(t, 1) if t else '0', 11, MUTED, 'text-anchor="middle"')]
    parts += [line(TX(0), top, TX(0), bottom, '#5B6B7C', 1.2),
              text((left + right) / 2, bottom + 36, 'θ, program units (secretory × ENDO interaction; negative = attenuation)', 12, INK, 'text-anchor="middle"')]
    for i, c in enumerate(COHORTS):
        y = top + 45 + i * 80; t = theta[c]
        x_shade = TX(max(tmin, -t['mde'] - 0.6))
        parts += [rect(x_shade, y - 22, TX(-t['mde']) - x_shade, 44, TINT[c]),
                  line(TX(-t['mde']), y - 22, TX(-t['mde']), y + 22, COLOR[c], 1.2, 'stroke-dasharray="3 3"'),
                  line(TX(t['ci_low']), y, TX(t['ci_high']), y, COLOR[c], 2.5),
                  line(TX(t['ci_low']), y - 7, TX(t['ci_low']), y + 7, COLOR[c], 2), line(TX(t['ci_high']), y - 7, TX(t['ci_high']), y + 7, COLOR[c], 2),
                  circle(TX(t['estimate']), y, 6, COLOR[c]),
                  text(cx0, y - 30, f"{NAME[c]}: θ {signed(t['estimate'])} [{signed(t['ci_low'])}, {signed(t['ci_high'])}], p {t['p']:.2f} · {t['secretory_endo']} ENDO / {t['secretory_non_endo']} non-ENDO secretory", 10.5, INK),
                  text(TX(-t['mde']) - 4, y + 36, f"detectable attenuation −{t['mde']:.2f}", 10.5, COLOR[c], 'text-anchor="end"')]
    parts += [text(cx0, 428, f"Gate AUC {theta['discovery']['auc']:.3f} / {theta['heca']['auc']:.3f}; secretory shift {theta['discovery']['shift']:.2f} / {theta['heca']['shift']:.2f} units.", 11.5, MUTED),
              text(cx0, 448, 'Shaded: detectable with 80% power. Label: no_detectable_difference.', 11.5, MUTED),
              text(cx0, LINES[0], 'Against each cohort’s own secretory shift, the intervals exclude:', 12.5),
              text(cx0, LINES[1], f"attenuation above {theta['discovery']['max_attenuation_share'] * 100:.1f}% (discovery) and {theta['heca']['max_attenuation_share'] * 100:.1f}% (HECA);", 11.5, MUTED),
              text(cx0, LINES[2], f"enhancement above {theta['discovery']['max_enhancement_share'] * 100:.1f}% and {theta['heca']['max_enhancement_share'] * 100:.1f}%.", 11.5, MUTED),
              text(cx0, LINES[3], 'Reading: no detectable program-level difference; the HECA bound is wide.', 11.5, INK)]

    parts += [line(35, 590, W - 35, 590),
              text(35, 610, 'Sources: reports/mde_analysis (v1 and phase MDE summaries, subsampling curves), reports/secretory_program (summaries), receipts/phase_positive_control_evaluation.json, reports/replication_v2_dossier_pool.tsv.', 11, MUTED),
              text(35, 628, 'Rules for every analysis were frozen and hashed before data; v1 artifacts unchanged. Receipt: receipts/summary_figure_follow_on_v1.json. Nothing here is disease evidence or a target claim.', 11, MUTED),
              '</g></svg>']
    return '\n'.join(parts) + '\n'


def main(root=ROOT):
    d = load_and_check(root)
    out = root / OUTPUT; out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(draw(d), encoding='utf-8')
    receipt = {'schema_version': 1, 'created_utc': datetime.now(timezone.utc).isoformat(),
               'purpose': 'Render one summary figure of the follow-on analyses from published summaries and receipts; no scientific model rerun',
               'input_sha256': {rel: sha(root / rel) for rel in INPUTS},
               'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'output': OUTPUT, 'output_sha256': sha(out),
               'numbers_drawn': {'disease_genes_q_lt_0_05': d['disease_sig'], 'phase_genes_q_lt_0_05': d['phase_sig'], 'v1_family_size': d['family'],
                                 'v1_B0_top20_transfer': d['v1_transfer'],
                                 'phase_B0_top20_transfer': {k: d['P2'][k] for k in ('same_sign', 'evaluable', 'selected', 'same_sign_and_external_q_lt_0_05')},
                                 'median_mde_fold': d['mde'], 'recall_at_11_per_arm': d['at11'], 'subsampling_curves_primary': d['curves'], 'theta': d['theta']},
               'new_scientific_analysis': False}
    (root / RECEIPT).parent.mkdir(parents=True, exist_ok=True)
    (root / RECEIPT).write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    return receipt


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', help='project root to read from and write into (default: this repository)')
    args = ap.parse_args()
    r = main(Path(args.root).resolve() if args.root else ROOT)
    print(json.dumps({k: r[k] for k in ('output', 'output_sha256', 'script_sha256')}))
