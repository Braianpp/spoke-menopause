# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Replot the saved Q2 cumulative DWPC values using the Figure 3C colors."""
import os
import csv
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/q3-matplotlib-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
STAGE = ROOT / 'audit_2026-09-30/excel_stage'
parser = argparse.ArgumentParser()
parser.add_argument('--layout', choices=['standalone', 'aligned'], default='standalone')
args = parser.parse_args()
aligned = args.layout == 'aligned'
stem = 'Figure3D_aligned_to_3C' if aligned else 'Figure3D_consistent_colors'

def read_csv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))
summary_path = STAGE / 'Q2_Target_Summary.csv'
metapath_path = STAGE / 'Q2_Metapaths_PM_Neuro.csv'
summary = {row['target']: row for row in read_csv(summary_path)}
metapaths = read_csv(metapath_path)
colors = json.loads((OUT / 'Figure3C_all27_validation.json').read_text())['target_colors']
order = ['Vascular dementia', 'Frontotemporal dementia', "Parkinson's disease"]
labels = ['Vascular dementia', 'Frontotemporal\ndementia', 'Parkinson’s disease']
records = []
for target in order:
    rows = [row for row in metapaths if row['target_name'] == target]
    total = math.fsum((float(row['dwpc']) for row in rows))
    saved = float(summary[target]['cumulative_DWPC'])
    assert math.isclose(total, saved, rel_tol=0, abs_tol=1e-12), (target, total, saved)
    records.append({'target': target, 'cumulative_DWPC': saved, 'sum_from_metapaths': total, 'color': colors[target]})
ink = '#162E3E'
plt.rcParams.update({'font.family': 'Arial', 'text.color': ink, 'axes.labelcolor': ink, 'xtick.color': ink, 'ytick.color': ink, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
if aligned:
    fig, ax = plt.subplots(figsize=(9.6, 8.6))
    fig.subplots_adjust(left=0.07, right=0.95, bottom=1.4 / 8.6, top=7.7 / 8.6)
    row_y = [6.95, 4.85, 2.25]
    bar_height = 0.92
else:
    fig, ax = plt.subplots(figsize=(10.5, 5.7))
    fig.subplots_adjust(left=0.345, right=0.97, bottom=0.18, top=0.85)
    row_y = list(range(3))
    bar_height = 0.55
bars = ax.barh(row_y, [row['cumulative_DWPC'] for row in records], height=bar_height, color=[row['color'] for row in records], edgecolor='#465761', linewidth=1.1)
ax.set_yticks([] if aligned else row_y, [] if aligned else labels, fontsize=22)
ax.tick_params(axis='y', length=0, pad=13)
ax.tick_params(axis='x', labelsize=24 if aligned else 20, pad=7)
ax.set_ylim((1.4, 7.7) if aligned else (2.58, -0.58))
ax.set_xlim(0, 0.06)
ax.set_xticks([0, 0.02, 0.04, 0.06], ['0.00', '0.02', '0.04', '0.06'])
ax.set_xlabel('Cumulative DWPC', fontsize=26 if aligned else 23, fontweight='bold', labelpad=14)
if not aligned:
    ax.set_title('Network connectivity strength', fontsize=23, fontweight='bold', loc='left', pad=20)
for side in ['top', 'right']:
    ax.spines[side].set_visible(False)
for side in ['left', 'bottom']:
    ax.spines[side].set_color('#71858C')
    ax.spines[side].set_linewidth(0.9)
for bar, row in zip(bars, records):
    ax.text(bar.get_width() + 0.0012, bar.get_y() + bar.get_height() / 2, f"{row['cumulative_DWPC']:.3f}", va='center', fontsize=24 if aligned else 20)
    assert to_hex(bar.get_facecolor()).upper() == row['color'].upper()
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
texts = [*ax.get_yticklabels(), *ax.get_xticklabels(), ax.xaxis.label, ax.title, ax._left_title, *ax.texts]
for text in texts:
    if not text.get_text():
        continue
    box = text.get_window_extent(renderer)
    assert box.x0 >= 0 and box.y0 >= 0 and (box.x1 <= fig.bbox.width) and (box.y1 <= fig.bbox.height), text.get_text()
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(OUT / f'{stem}.{ext}', dpi=600, facecolor='white')
fig.savefig(OUT / ('Figure3D_aligned_preview.png' if aligned else 'Figure3D_preview.png'), dpi=100, facecolor='white')
fig.savefig(OUT / ('Figure3D_aligned_panel_size.png' if aligned else 'Figure3D_panel_size.png'), dpi=700 / fig.get_figwidth(), facecolor='white')
row_centers = [float(fig.transFigure.inverted().transform(ax.transData.transform((0, y)))[1]) for y in row_y]
if aligned:
    assert all((math.isclose(actual, expected / 8.6, abs_tol=1e-12) for actual, expected in zip(row_centers, [6.95, 4.85, 2.25])))
plt.close(fig)
validation = {'records': records, 'all_colors_match_Figure3C': True, 'all_values_match_saved_summary_and_metapath_sums': True, 'all_text_within_canvas': True, 'disease_order_top_to_bottom': order, 'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [summary_path, metapath_path]}, 'original_plot_script': '${SPOKE_LOCAL_INPUT}', 'changes': 'Disease colors match Figure 3C at full opacity; disease labels visible; Arial text and vector exports.', 'visual_review': 'Pending'}
if aligned:
    validation.update({'layout': 'Companion to Figure 3C; equal uncropped image height, top aligned.', 'canvas_inches': [9.6, 8.6], 'row_centers_fraction_from_bottom': row_centers, 'row_centers_match_Figure3C': True, 'changes': 'Increased vertical spacing to match all three Figure 3C target-node positions; no repeated disease labels or title in the assembly version; exact disease colors and DWPC values preserved.'})
(OUT / ('Figure3D_aligned_validation.json' if aligned else 'Figure3D_validation.json')).write_text(json.dumps(validation, indent=2) + '\n')
(OUT / ('Figure3D_aligned_caption.txt' if aligned else 'Figure3D_caption.txt')).write_text('d, Cumulative degree-weighted path count (DWPC), summed across evaluated metapaths from premature menopause to each target disease. Disease colors match panel c: vascular dementia, purple; frontotemporal dementia, blue; Parkinson’s disease, green. Values beside the bars are rounded to three decimal places; bar lengths use the saved full-precision values (Supplementary Table 3, sections 3a–b).\n')
print(json.dumps({'outputs': stem + '.{png,pdf,svg}', 'records': records}, indent=2))
