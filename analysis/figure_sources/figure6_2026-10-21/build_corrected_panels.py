# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Replot Figure 6c/d from the verified submission tables, without new queries."""
import os
from pathlib import Path
from collections import defaultdict
from zipfile import ZipFile
import hashlib
import json
import os
import xml.etree.ElementTree as ET
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/spoke-figure6-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
import numpy as np
OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent
WORKBOOK = ROOT / 'outputs/supplementary_submission_2026-10-08/Supplementary_Tables_Submission.xlsx'
payload = json.loads((ROOT / 'supplementary_tables_2026-10-01/payload_submission.json').read_text())
tables = {int(t['number']): t for t in payload['tables']}
NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
REL = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'
with ZipFile(WORKBOOK) as z:
    strings = []
    if 'xl/sharedStrings.xml' in z.namelist():
        strings = [''.join((t.text or '' for t in si.iter(NS + 't'))) for si in ET.fromstring(z.read('xl/sharedStrings.xml'))]
    rels = {r.get('Id'): r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    sheet_rows = {}
    for sh in ET.fromstring(z.read('xl/workbook.xml')).find(NS + 'sheets'):
        if not sh.get('name').startswith(('S32_', 'S37_', 'S39_', 'S40_')):
            continue
        target = rels[sh.get(REL + 'id')]
        target = target.lstrip('/') if target.startswith('/') else 'xl/' + target
        rows = []
        for row in ET.fromstring(z.read(target)).iter(NS + 'row'):
            values = {}
            for cell in row:
                col = ''.join((x for x in cell.get('r') if x.isalpha()))
                ix = 0
                for x in col:
                    ix = ix * 26 + ord(x) - 64
                v = cell.find(NS + 'v')
                if cell.get('t') == 's':
                    val = strings[int(v.text)]
                elif cell.get('t') == 'inlineStr':
                    val = ''.join((x.text or '' for x in cell.iter(NS + 't')))
                elif cell.get('t') in ('str', 'e'):
                    val = v.text if v is not None else None
                elif v is not None:
                    val = float(v.text)
                else:
                    val = None
                values[ix - 1] = val
            rows.append([values.get(i) for i in range(max(values, default=-1) + 1)])
        sheet_rows[int(sh.get('name')[1:3])] = rows

def same_row(a, b):
    if len(b) < len(a):
        b = b + [None] * (len(a) - len(b))
    for x, y in zip(a, b):
        if x is None and y in (None, ''):
            continue
        if isinstance(x, (int, float)) and isinstance(y, (int, float)):
            if not np.isclose(x, y, atol=1e-14, rtol=1e-12):
                return False
        elif x != y:
            return False
    return True
for n in (32, 37, 39, 40):
    blocks = tables[n]['blocks'][:1] if n == 37 else tables[n]['blocks']
    for block in blocks:
        for row in block['rows']:
            assert any((same_row(row, candidate) for candidate in sheet_rows[n])), (n, row)
route_rows = tables[37]['blocks'][0]['rows']
compound_ppr = tables[32]['blocks'][0]['rows']
histogram = tables[32]['blocks'][1]['rows']
assert len(histogram) == 50
assert sum((row[2] for row in histogram)) == 550596
route_data = {(r[0], r[2]): r[3] for r in route_rows if r[1] == "Alzheimer's disease"}
categories = ['Gene regulation', 'Other', 'Protein interaction', 'Genetic variant association', 'Biomarker overlap', 'Therapeutic overlap', 'Disease hierarchy']
labels = ['Gene\nregulation', 'Other', 'Protein\ninteraction', 'Genetic variant\nassociation', 'Biomarker\noverlap', 'Therapeutic\noverlap', 'Disease\nhierarchy']
genes = ['ESR1', 'CYP2B6', 'CYP2C19', 'CYP2C9', 'SULT1A1', 'UGT1A1', 'TP53', 'NFE2L2']
compounds = ['Estradiol', 'Estradiol acetate', 'Ethinyl estradiol', 'Equilin', 'Estrone', 'Bazedoxifene']
edges = tables[40]['blocks'][0]['rows']
connections = defaultdict(set)
for row in edges:
    connections[row[3], row[1].lower()].add(row[4])
code_for = {frozenset(): 0, frozenset(['BINDS']): 1, frozenset(['INTERACTS']): 2, frozenset(['BINDS', 'INTERACTS']): 3, frozenset(['DOWNREGULATES']): 5, frozenset(['INTERACTS', 'UPREGULATES', 'DOWNREGULATES']): 6}
matrix = np.array([[code_for[frozenset(connections[g, c.lower()])] for c in compounds] for g in genes])
saved = tables[39]['blocks'][1]
assert saved['headers'][1:] == compounds
assert [r[0] for r in saved['rows']] == genes
assert np.array_equal(matrix, np.array([r[1:] for r in saved['rows']]))
assert connections['TP53', 'estradiol'] == {'INTERACTS', 'UPREGULATES', 'DOWNREGULATES'}
assert np.count_nonzero(matrix) == 19
assert len(edges) == 28
INK = '#162E3E'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 14, 'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none', 'pdf.fonttype': 42})

def panel_b(fig, box):
    x0, y0, w, h = box
    ax = fig.add_axes([x0 + 0.15 * w, y0 + 0.31 * h, 0.82 * w, 0.55 * h])
    ax.bar([r[0] for r in histogram], [r[2] for r in histogram], width=[r[1] - r[0] for r in histogram], align='edge', color='#D9D9D9', edgecolor='white', linewidth=0.6, zorder=2)
    colors = {'Tamoxifen': '#1F77B4', 'Leuprolide': '#FF7F0E', 'Letrozole': '#2CA02C', 'Goserelin': '#D62728'}
    for row in compound_ppr:
        ax.axvline(np.log10(row[2] + 1e-10), color=colors[row[0]], ls=(0, (5, 3)), lw=1.9, label=row[0], zorder=3)
    ax.set_xticks(np.arange(-9, -2, 1))
    ax.set_xlim(-9, histogram[-1][1])
    ax.set_ylim(0, 65000)
    ax.set_yticks([0, 20000, 40000, 60000], ['0', '20,000', '40,000', '60,000'])
    ax.set_xlabel('$\\log_{10}$(PageRank + $10^{-10}$)', fontsize=16, fontweight='bold', labelpad=8)
    ax.set_ylabel('Compound count', fontsize=16, fontweight='bold', labelpad=8)
    ax.grid(axis='y', color='#E6EBEF', linewidth=0.7, zorder=0)
    for s in ['left', 'bottom']:
        ax.spines[s].set_color('#9AA9B1')
    fig.text(x0 + 0.035 * w, y0 + 0.975 * h, 'B. Q8 Antiestrogen therapies–Dementia', fontsize=17, fontweight='bold', ha='left', va='top')
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(x0 + 0.53 * w, y0 + 0.035 * h), bbox_transform=fig.transFigure, ncol=2, frameon=False, fontsize=14, columnspacing=2)
    return ax

def panel_c(fig, box):
    x0, y0, w, h = box
    ax = fig.add_axes([x0 + 0.13 * w, y0 + 0.27 * h, 0.84 * w, 0.59 * h])
    x = np.arange(len(categories))
    bw = 0.35
    for delta, compound, color in [(-bw / 2, 'Estradiol', '#F4A261'), (bw / 2, 'Estrone', '#E63946')]:
        ax.bar(x + delta, [route_data[compound, c] for c in categories], bw, label=compound, color=color, edgecolor='white', linewidth=0.4, zorder=3)
    ax.set_xticks(x, labels, fontsize=12.5, rotation=32, ha='right', rotation_mode='anchor')
    ax.set_ylim(0, 0.4)
    ax.set_yticks(np.arange(0, 0.401, 0.1))
    ax.set_ylabel('Cumulative DWPC', fontsize=16, fontweight='bold', labelpad=7)
    ax.grid(axis='y', color='#E6EBEF', linewidth=0.7, zorder=0)
    ax.legend(loc='upper right', frameon=False, fontsize=14)
    for s in ['left', 'bottom']:
        ax.spines[s].set_color('#9AA9B1')
    fig.text(x0 + 0.035 * w, y0 + 0.975 * h, 'C. Q9 Estrogen compounds–Alzheimer’s disease', fontsize=17, fontweight='bold', ha='left', va='top')
    return ax

def panel_d(fig, box):
    x0, y0, w, h = box
    ax = fig.add_axes([x0 + 0.18 * w, y0 + 0.39 * h, 0.79 * w, 0.47 * h])
    colors = ['#F5F7F9', '#3B82F6', '#10B981', '#8B5CF6', '#F97316', '#EF4444', '#1E3A5F']
    ax.imshow(matrix, cmap=ListedColormap(colors), norm=BoundaryNorm(np.arange(-0.5, 7.5, 1), 7), aspect='auto')
    ax.set_yticks(range(8), genes, fontsize=14, fontstyle='italic')
    ax.set_xticks(range(6), ['Estradiol', 'Estradiol\nacetate', 'Ethinyl\nestradiol', 'Equilin', 'Estrone', 'Bazedoxifene'], fontsize=13, rotation=38, ha='right', rotation_mode='anchor')
    ax.set_xticks(np.arange(-0.5, 6, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 8, 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=1.4)
    ax.tick_params(which='both', length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.text(x0 + 0.035 * w, y0 + 0.975 * h, 'D. Q10 Pharmacogenetics of estrogen compounds', fontsize=17, fontweight='bold', ha='left', va='top')
    entries = [(1, 'Binds protein'), (2, 'Interacts with protein'), (3, 'Binding + protein interaction'), (5, 'Downregulates gene'), (6, 'Protein interaction +\nup- and downregulation'), (0, 'No retrieved relationship')]
    handles = [Patch(facecolor=colors[k], edgecolor='#AAB5BE', linewidth=0.5, label=label) for k, label in entries]
    fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(x0 + 0.025 * w, y0 + 0.01 * h), bbox_transform=fig.transFigure, ncol=2, frameon=False, fontsize=12.5, handlelength=1.1, columnspacing=1.1, labelspacing=0.65)
    return ax
figures = []
for name, draw in [('Figure6B_corrected', panel_b), ('Figure6C_corrected', panel_c), ('Figure6D_corrected', panel_d)]:
    fig = plt.figure(figsize=(7.6, 5.8), facecolor='white')
    draw(fig, (0, 0, 1, 1))
    fig.savefig(OUT / (name + '.png'), dpi=600, facecolor='white')
    fig.savefig(OUT / (name + '.svg'), facecolor='white')
    fig.savefig('/private/tmp/' + name + '_preview.png', dpi=130, facecolor='white')
    figures.append(fig)
fig = plt.figure(figsize=(15.2, 5.8), facecolor='white')
panel_c(fig, (0, 0, 0.5, 1))
panel_d(fig, (0.5, 0, 0.5, 1))
fig.savefig(OUT / 'Figure6CD_corrected.png', dpi=600, facecolor='white')
fig.savefig(OUT / 'Figure6CD_corrected.svg', facecolor='white')
fig.savefig('/private/tmp/Figure6CD_corrected_preview.png', dpi=100, facecolor='white')
for f in figures + [fig]:
    f.canvas.draw()
    renderer = f.canvas.get_renderer()
    fw, fh = (f.bbox.width, f.bbox.height)
    for txt in f.findobj(matplotlib.text.Text):
        if not txt.get_visible() or not txt.get_text():
            continue
        bb = txt.get_window_extent(renderer)
        assert bb.x0 >= -1 and bb.y0 >= -1 and (bb.x1 <= fw + 1) and (bb.y1 <= fh + 1), txt.get_text()
    plt.close(f)
report = {'source': str(WORKBOOK), 'source_sha256': hashlib.sha256(WORKBOOK.read_bytes()).hexdigest(), 'tables': [32, 37, 39, 40], 'input_rows_match_submission_workbook': True, 'bar_values_unchanged': True, 'matrix_cells_unchanged': True, 'histogram_bins_unchanged': 50, 'histogram_count': 550596, 'panel_b_full_saved_histogram_range': True, 'panel_b_original_four_ppr_values_unchanged': True, 'edge_records': len(edges), 'connected_cells': int(np.count_nonzero(matrix)), 'tp53_estradiol_relationships': sorted(connections['TP53', 'estradiol']), 'target_in_panel_c': "Alzheimer's disease", 'figure_text_within_bounds': True}
(OUT / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
