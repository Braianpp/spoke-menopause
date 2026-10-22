# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Rebuild Figure 3 with shared genes and separately analyzed PD mitochondria."""
import os
import csv
import hashlib
import json
import math
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/q2-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Wedge, FancyBboxPatch, PathPatch, Patch, Rectangle
from matplotlib.path import Path as MP
from matplotlib.transforms import Bbox
OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
OLD = ROOT / 'figures/figure3_2026-10-01'
PD = ROOT / 'analyses/q2_pd_mitochondria_2026-10-06'

def read(p):
    with p.open(newline='') as f:
        return list(csv.DictReader(f))
meta = json.loads((OLD / 'Figure3C_all27_validation.json').read_text())
shared_terms = read(OLD / 'Figure3B_colored_terms_and_genes.csv')
pd_terms = read(PD / 'PD_165_selected_mitochondrial_terms.csv')
assert len(shared_terms) == 20 and len(pd_terms) == 3
autophagy = next((r for r in pd_terms if 'GO:1903146' in r['term']))
pd_genes = ['PINK1', 'PRKN', 'PARK7', 'FBXO7', 'HTRA2', 'HTT']
assert set(pd_genes) == set(autophagy['genes'].split(';'))
all_genes = {r['gene_name']: r for r in read(ROOT / 'audit_2026-09-30/excel_stage/Q2_IntermediateGenes_242.csv')}
assert all((all_genes[g]['targets'] == "Parkinson's disease" for g in pd_genes))
assert not set(pd_genes) & set(meta['gene_positions'])
families = {'Glial development': '#96BDB7', 'mtDNA maintenance': '#598EAF', 'Superoxide regulation': '#C59B55', 'Microglial activation': '#B77798', 'Mitochondrial quality control': '#287F8E'}
ink, grey, edge = ('#162E3E', '#BAC3CA', '#A1AEB5')
colors = meta['target_colors']
plt.rcParams.update({'font.family': 'Arial', 'text.color': ink, 'axes.labelcolor': ink, 'xtick.color': ink, 'ytick.color': ink, 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold', 'mathtext.default': 'regular'})
W, H = (22.8, 26.25)
fig = plt.figure(figsize=(W, H))
canvas = fig.add_axes([0, 0, 1, 1])
canvas.set(xlim=(0, W), ylim=(0, H))
canvas.axis('off')
texts = []

def title(x, y, s):
    texts.append(canvas.text(x, y, s, fontsize=30, color='black', va='center'))

def label(x, y, s, fs=22, **kw):
    a = canvas.text(x, y, s, fontsize=fs, va='center', **kw)
    texts.append(a)
    return a
title(0.2, 24.95, 'A. Network proximity')
title(13.45, 24.95, 'B. Biological processes')
label(0.65, 24.35, 'Premature menopause to diseases in SPOKE', 23)
label(13.45, 24.35, 'Shared intermediates (27 genes)', 23, fontweight='bold')
proximity = read(ROOT / 'audit_2026-09-30/excel_stage/Q2_DiseaseProximity_2422.csv')
by_id = {r['disease_id']: r for r in proximity}
values = np.array([int(r['shared_3hop']) for r in proximity])
assert len(values) == 2422
bins = np.arange(0, values.max() + 6, 5)
counts, edges = np.histogram(values, bins)
assert counts.sum() == 2422
a = fig.add_axes([1.25 / W, 17.05 / H, 10.5 / W, 6.65 / H])
a.hist(values, bins=bins, orientation='horizontal', color='#D9D9D9', edgecolor='white', lw=0.6, zorder=2)
a.set(xlim=(0, 1500), ylim=(0, 70), xticks=[0, 500, 1000, 1500], yticks=[0, 20, 40, 60, 70])
a.tick_params(labelsize=23, pad=7)
a.set_xlabel('Number of diseases', fontsize=26, fontweight='bold', labelpad=13)
a.set_ylabel('Shared 3-hop gene intermediates', fontsize=25, fontweight='bold', labelpad=18)
a.grid(axis='x', color='#E5E9EC', lw=0.8)
a.set_axisbelow(True)
a.text(0.04, 0.96, 'Detail: 0–70 intermediates', transform=a.transAxes, fontsize=20, va='top')
handles = [Patch(facecolor='#D9D9D9', edgecolor='white', label='All diseases')]
markers = []
for did, name in [('DOID:8725', 'Vascular dementia'), ('DOID:9255', 'Frontotemporal dementia'), ('DOID:14330', "Parkinson's disease"), ('DOID:10652', "Alzheimer's disease")]:
    v = int(by_id[did]['shared_3hop'])
    color = colors.get(name, '#E67E22')
    a.axhline(v, color=color, lw=3, ls='--', zorder=4)
    handles.append(Line2D([0], [0], color=color, lw=3, ls='--', label=name.replace("'", '’')))
    markers.append({'disease': name, 'shared_3hop': v, 'color': color})
for did in ['DOID:332', 'DOID:12217', 'DOID:678', 'DOID:12858']:
    a.axhline(int(by_id[did]['shared_3hop']), color='#8C99A4', lw=1.7, zorder=3)
handles.append(Line2D([0], [0], color='#8C99A4', lw=1.7, label='Other neurodegenerative references'))
legend = fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(6.2 / W, 13.6 / H), frameon=False, fontsize=21, handlelength=1.8, handletextpad=0.7, labelspacing=0.5, borderaxespad=0)
texts.extend(legend.get_texts())
inset = fig.add_axes([1.25 / W, 14.03 / H, 4.25 / W, 1.75 / H])
nz = counts > 0
inset.barh((edges[:-1] + edges[1:])[nz] / 2, counts[nz] - 0.8, left=0.8, height=5, color='#ACB8C0', edgecolor='none')
inset.set(xscale='log', xlim=(0.8, 2000), ylim=(0, 1120))
inset.set_xticks([1, 10, 100, 1000], ['1', '10', '100', '1,000'])
inset.set_yticks([0, 500, 1000], ['0', '500', '1,000'])
inset.tick_params(labelsize=18, pad=4)
inset.set_title('Full range: all 2,422 diseases', fontsize=20, loc='left', pad=9)
inset.set_xlabel('Number of diseases (log scale)', fontsize=19, labelpad=6)
widths = np.array([1.85 if '\n' in r['plot_label'] else 1.0 for r in shared_terms])
centers = np.cumsum(widths) - widths / 2
extent = sum(widths)
b = fig.add_axes([18.05 / W, 15.55 / H, 4.05 / W, 8.15 / H])
vs = [-math.log10(float(r['Adjusted P-value'])) for r in shared_terms]
assert all((math.isclose(v, float(r['minus_log10_adjusted_P']), abs_tol=1e-12) for v, r in zip(vs, shared_terms)))
b.barh(centers, vs, height=0.59, color=[r['color'] for r in shared_terms], edgecolor='#71858C', lw=0.65)
b.set_yticks(centers, [r['plot_label'] for r in shared_terms], fontsize=22)
b.set(xlim=(0, 7.9), ylim=(extent + 0.1, -0.1), xticks=[0, 2, 4, 6])
b.tick_params(axis='y', length=0, pad=10)
b.tick_params(axis='x', labelsize=20, pad=6)
for y, v, r in zip(centers, vs, shared_terms):
    b.text(v + 0.09, y, 'n=' + r['n_overlap'], va='center', fontsize=18)
b.set_xlabel('$-\\log_{10}$(adjusted P value)', fontsize=22, fontweight='bold', labelpad=10)
label(13.45, 14.45, 'Parkinson’s intermediates (165 genes)', 23, fontweight='bold', color=families['Mitochondrial quality control'])
label(13.45, 14.05, 'Selected mitochondrial processes', 19)
bp = fig.add_axes([18.05 / W, 11.73 / H, 4.05 / W, 2.18 / H])
plabels = ['Regulation of mitochondrial\nautophagy', 'Regulation of\nmitochondrial fission', 'Regulation of mitochondrial\nmembrane potential']
pvs = [-math.log10(float(r['adjusted_p_value'])) for r in pd_terms]
bp.barh([0, 1, 2], pvs, height=0.53, color=families['Mitochondrial quality control'], edgecolor='#46646B', lw=0.65)
bp.set_yticks([0, 1, 2], plabels, fontsize=22)
bp.set(xlim=(0, 7.9), ylim=(2.55, -0.55), xticks=[0, 2, 4, 6])
bp.tick_params(axis='y', length=0, pad=10)
bp.tick_params(axis='x', labelsize=20, pad=6)
bp.set_xlabel('$-\\log_{10}$(adjusted P value)', fontsize=22, fontweight='bold', labelpad=10)
for y, v, r in zip(range(3), pvs, pd_terms):
    bp.text(v + 0.09, y, 'n=' + r['overlap_count'], va='center', fontsize=18)
for axis in [b, bp]:
    axis.set_axisbelow(True)
    axis.grid(axis='x', color='#E5E9EC', lw=0.7)
upper_shift = 0.65
for artist in texts:
    if artist.axes is canvas:
        x, y = artist.get_position()
        artist.set_position((x, y + upper_shift))
for axis in [a, inset, b, bp]:
    pos = axis.get_position()
    axis.set_position([pos.x0, pos.y0 + upper_shift / H, pos.width, pos.height])
legend.set_bbox_to_anchor((6.2 / W, (13.6 + upper_shift) / H))
title(0.2, 10.75, 'C. Shared and PD mitochondrial intermediates')
title(13.45, 10.75, 'D. Network connectivity strength')
cx, cy, cw, ch = (0.2, 0.2, 13.1, 10.1)
c = fig.add_axes([cx / W, cy / H, cw / W, ch / H])
c.set(xlim=(0, cw), ylim=(0, ch), aspect='equal')
c.axis('off')
ctexts = []
gene_labels = []
gene_markers = {}

def ct(x, y, s, fs=23, **kw):
    artist = c.text(x, y, s, fontsize=fs, color=ink, **kw)
    ctexts.append(artist)
    return artist

def curve(points):
    c.add_patch(PathPatch(MP(points, [MP.MOVETO] + [MP.CURVE4] * 3), facecolor='none', edgecolor=edge, lw=1.55, zorder=0))
source_y = 5.7
target_y = [8.35, 6.2, 2.85]
for end in [7.2, 4.05, 2.15]:
    curve([(2.12, source_y), (2.5, source_y), (2.3, end), (2.7, end)])
for start, end in [(8.1, 8.35), (6.9, 6.2), (4.35, 6.2), (3.7, 2.85), (2.15, 2.85)]:
    curve([(9.1, start), (9.44, start), (9.18, end), (9.54, end)])
c.add_patch(Circle((1.16, source_y), 0.97, facecolor=meta['source_color'], edgecolor='white', lw=1, zorder=2))
ctexts.append(c.text(1.16, source_y, 'Premature\nmenopause', ha='center', va='center', fontsize=23, color='white', zorder=3))
ct(1.16, 9.91, 'Source', 22, ha='center', va='center')
ct(5.78, 9.91, '27 shared + 6 selected PD genes', 23, ha='center', va='center', fontweight='bold')
ct(10.85, 9.91, 'Targets', 22, ha='center', va='center')
groups = [(meta['VaD_FTD_genes'], 5.1, 4.35, 'Vascular dementia + FTD · 21 genes'), (meta['FTD_PD_genes'], 3.3, 1.49, 'FTD + Parkinson’s disease · 6 genes'), (pd_genes, 1.4, 1.49, 'PD only · 6 mitochondrial autophagy genes')]
gene_records = []
for gi, (gs, y0, height, heading) in enumerate(groups):
    c.add_patch(FancyBboxPatch((2.7, y0), 6.4, height, boxstyle='round,pad=0.04,rounding_size=0.1', facecolor='#EDF6F8' if gi == 2 else '#F8FAFB', edgecolor=families['Mitochondrial quality control'] if gi == 2 else '#D4DEE3', linewidth=1.0, zorder=1))
    ct(2.9, y0 + height - 0.3, heading, 20.5, ha='left', va='center', fontweight='bold')
    for j, g in enumerate(gs):
        x = [2.95, 5.25, 7.25][j % 3]
        y = y0 + height - 0.72 - 0.54 * (j // 3)
        ann = ['Mitochondrial quality control'] if gi == 2 else meta['highlight_memberships'][g]
        cols = [co for fa, co in families.items() if fa in ann] or [grey]
        for k, co in enumerate(cols):
            c.add_patch(Wedge((x, y), 0.105, 90 + 360 * k / len(cols), 90 + 360 * (k + 1) / len(cols), facecolor=co, edgecolor='none', zorder=3))
        marker = Circle((x, y), 0.105, facecolor='none', edgecolor='#7E939F', lw=0.6, zorder=4)
        c.add_patch(marker)
        gene_markers[g] = marker
        gene_labels.append(ct(x + 0.22, y, g, 26, ha='left', va='center', fontstyle='italic', fontweight='bold' if gi == 2 or g in {'POLG', 'POLG2', 'TWNK'} else 'normal', zorder=4))
        gene_records.append({'gene': g, 'targets': all_genes[g]['targets'], 'display_group': heading, 'selection': 'All contributors to GO:1903146 in full-PD enrichment' if gi == 2 else 'All genes shared by two target sets', 'annotations': ';'.join(ann), 'x': x, 'y': y})
targets = ['Vascular dementia', 'Frontotemporal dementia', "Parkinson's disease"]
tlabels = ['Vascular\ndementia', 'Frontotemporal\ndementia', 'Parkinson’s\ndisease']
target_labels = []
for t, lab, y in zip(targets, tlabels, target_y):
    c.add_patch(Circle((9.7, y), 0.14, facecolor=colors[t], edgecolor='white', lw=0.8, zorder=3))
    target_labels.append(ct(9.98, y, lab, 24, ha='left', va='center'))
for x, y, lab in [(0.55, 0.84, 'Glial development'), (4.45, 0.84, 'Superoxide regulation'), (8.67, 0.84, 'Other shared genes'), (0.55, 0.34, 'mtDNA maintenance'), (4.45, 0.34, 'Microglial activation'), (8.67, 0.34, 'Mitochondrial quality control')]:
    c.add_patch(Rectangle((x, y - 0.1), 0.2, 0.2, facecolor=families.get(lab, grey), edgecolor='#7E939F', lw=0.6))
    ct(x + 0.34, y, lab, 19.5, ha='left', va='center')
d = fig.add_axes([14.1 / W, (cy + 1.9) / H, 8.0 / W, 7.2 / H])
d.set(xlim=(0, 0.06), ylim=(1.9, 9.1))
d.set_yticks([])
d.set_xticks([0, 0.02, 0.04, 0.06], ['0.00', '0.02', '0.04', '0.06'])
d.tick_params(labelsize=24, pad=7)
d.set_xlabel('Cumulative DWPC', fontsize=26, fontweight='bold', labelpad=14)
summary = {r['target']: r for r in read(ROOT / 'audit_2026-09-30/excel_stage/Q2_Target_Summary.csv')}
metapaths = read(ROOT / 'audit_2026-09-30/excel_stage/Q2_Metapaths_PM_Neuro.csv')
scores = []
for t, y in zip(targets, target_y):
    value = float(summary[t]['cumulative_DWPC'])
    assert math.isclose(value, math.fsum((float(r['dwpc']) for r in metapaths if r['target_name'] == t)), abs_tol=1e-12)
    d.barh(y, value, height=0.82, color=colors[t], edgecolor='#465761', lw=1.1)
    d.text(value + 0.0012, y, f'{value:.3f}', fontsize=24, va='center')
    scores.append({'target': t, 'DWPC': value})
for axis in [a, inset, b, bp, d]:
    for side in ['top', 'right']:
        axis.spines[side].set_visible(False)
    for side in ['bottom', 'left']:
        axis.spines[side].set_color('#71858C')
        axis.spines[side].set_linewidth(0.8)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
check_texts = texts + ctexts
for axis in [a, inset, b, bp, d]:
    check_texts += [*axis.get_xticklabels(), *axis.get_yticklabels(), axis.xaxis.label, axis.yaxis.label, *axis.texts]
check_texts.append(inset._left_title)
for t in check_texts:
    bb = t.get_window_extent(renderer)
    assert bb.x0 >= 0 and bb.y0 >= 0 and (bb.x1 <= fig.bbox.width) and (bb.y1 <= fig.bbox.height), ('Outside figure', t.get_text(), bb)
for axis in [b, bp]:
    boxes = [t.get_window_extent(renderer) for t in axis.get_yticklabels()]
    assert all((boxes[i].y0 > boxes[i + 1].y1 for i in range(len(boxes) - 1))), 'B labels overlap'
    assert min((bb.x0 for bb in boxes)) > a.bbox.x1 + 10, 'B encroaches on A'
for i, t in enumerate(gene_labels):
    bb = t.get_window_extent(renderer)
    for u in gene_labels[i + 1:]:
        assert not bb.overlaps(u.get_window_extent(renderer)), (t.get_text(), u.get_text())
    for g, m in gene_markers.items():
        if g != t.get_text():
            assert not bb.overlaps(m.get_window_extent(renderer)), (t.get_text(), g)
for t in target_labels:
    assert t.get_window_extent(renderer).x1 < d.bbox.x0 - 20, 'C target touches D'
for y in target_y:
    assert math.isclose(c.transData.transform((9.7, y))[1], d.transData.transform((0, y))[1], abs_tol=1e-09)
assert bp.xaxis.label.get_window_extent(renderer).y0 > texts[-1].get_window_extent(renderer).y1 + 8, 'B axis label touches D title'
assert len(gene_records) == 33 and len({r['gene'] for r in gene_records}) == 33
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(OUT / f'Figure3_revised.{ext}', dpi=300, facecolor='white')
fig.savefig(OUT / 'Figure3_preview.png', dpi=90, facecolor='white')
for name, bounds in [('B', (12.15, 11.5, 22.8, H)), ('C', (0, 0, 13.65, 11.12)), ('CD', (0, 0, 22.8, 11.12))]:
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(OUT / f'Figure3{name}_revised.{ext}', dpi=300, facecolor='white', bbox_inches=Bbox.from_extents(*bounds), pad_inches=0)
plt.close(fig)
with (OUT / 'Figure3C_displayed_genes.csv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(gene_records[0]))
    writer.writeheader()
    writer.writerows(gene_records)
validation = {'genes_shown': 33, 'shared_genes_retained': 27, 'PD_autophagy_genes_added': pd_genes, 'selected_PD_term': 'GO:1903146', 'PD_enrichment_input_genes': 165, 'shared_enrichment_input_genes': 27, 'shared_terms_retained': 20, 'PD_terms_added': 3, 'gene_label_and_marker_overlap_check': 'passed', 'all_text_in_bounds': True, 'disease_rows_aligned_with_D': True, 'connectors': 'Summarized membership in extracted paths, not individual SPOKE edges', 'A_selected_markers': markers, 'D_scores': scores, 'source_hashes': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [OLD / 'Figure3B_colored_terms_and_genes.csv', PD / 'PD_165_selected_mitochondrial_terms.csv', ROOT / 'audit_2026-09-30/excel_stage/Q2_IntermediateGenes_242.csv']}, 'visual_review': 'Pending', 'font': 'Arial', 'canvas_inches': [W, H]}
(OUT / 'Figure3_validation.json').write_text(json.dumps(validation, indent=2) + '\n')
print(json.dumps({'figure': 'Figure3_revised', 'genes_shown': 33, 'PD_input': 165, 'shared_input': 27, 'geometry_checks': 'passed'}))
