# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Plot the saved representative Q3 enrichment terms; performs no web requests."""
import os
import argparse
import csv
import json
import math
import os
from pathlib import Path
ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--layout', choices=['compact', 'tall'], default='compact')
args = parser.parse_args()
tall = args.layout == 'tall'
OUT = ROOT.parent / 'figure4c_2026-10-01' if tall else ROOT
OUT.mkdir(exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/q3-matplotlib-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
rows = list(csv.DictReader((ROOT / 'Figure4C_plotted_terms.csv').open()))
plt.rcParams.update({'font.family': 'Arial', 'font.size': 16, 'text.color': '#162E3E', 'axes.labelcolor': '#162E3E', 'xtick.color': '#162E3E', 'ytick.color': '#162E3E', 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
colors = {'GO_Biological_Process_2023': '#96BDB7', 'KEGG_2021_Human': '#8E7DA7', 'Reactome_2022': '#69A3C0'}
fig, ax = plt.subplots(figsize=(12, 13 if tall else 10))
fig.subplots_adjust(left=0.56, right=0.965, top=0.975 if tall else 0.9, bottom=0.08 if tall else 0.145)
y = list(range(len(rows)))
values = [-math.log10(float(r['Adjusted P-value'])) for r in rows]
bars = ax.barh(y, values, height=0.52 if tall else 0.58, color=[colors[r['Library']] for r in rows], edgecolor='#516F77', linewidth=0.7)
ax.set_yticks(y, labels=[r['Plot_label'] for r in rows], fontsize=22 if tall else 18)
ax.tick_params(axis='y', length=0, pad=11)
ax.tick_params(axis='x', labelsize=17 if tall else 16, length=4)
ax.set_ylim(len(rows) - 0.4, -0.65)
ax.set_xlim(0, 15)
ax.set_xticks([0, 3, 6, 9, 12, 15])
ax.set_xlabel('$-\\log_{10}$(adjusted P value)', fontsize=19 if tall else 18, fontweight='bold' if tall else 'normal', labelpad=12)
ax.set_axisbelow(True)
ax.grid(axis='x', color='#DCE3E6', linewidth=0.7)
for spine in ['left', 'bottom']:
    ax.spines[spine].set_color('#9EAEB4')
    ax.spines[spine].set_linewidth(0.8)
for bar, r in zip(bars, rows):
    ax.text(bar.get_width() + 0.2, bar.get_y() + bar.get_height() / 2, 'n=' + r['Overlap_count'], va='center', fontsize=16.5 if tall else 15)
legend_entries = [('KEGG_2021_Human', 'KEGG'), ('GO_Biological_Process_2023', 'GO Biological Process'), ('Reactome_2022', 'Reactome')]
handles = [Patch(facecolor=colors[k], edgecolor='#516F77', label=label) for k, label in legend_entries]
if tall:
    ax.legend(handles=handles, loc='lower right', bbox_to_anchor=(1.035, 0.012), frameon=False, fontsize=17, handlelength=1.2, borderpad=0.3, labelspacing=0.6)
else:
    fig.text(0.5, 0.965, 'Enrichment of the 15 symptom-linked genes', ha='center', va='top', fontsize=22, fontweight='bold')
    fig.text(0.5, 0.925, 'Selected significant terms', ha='center', va='top', fontsize=17)
    fig.legend(handles=[handles[1], handles[0], handles[2]], loc='lower center', bbox_to_anchor=(0.5, 0.025), ncol=3, frameon=False, fontsize=16, handlelength=1.2, columnspacing=1.6)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
texts = [*ax.get_yticklabels(), *ax.get_xticklabels(), ax.xaxis.label, *ax.texts]
for text in texts:
    box = text.get_window_extent(renderer)
    assert box.x0 >= 0 and box.y0 >= 0 and (box.x1 <= fig.bbox.width) and (box.y1 <= fig.bbox.height), text.get_text()
for ext, dpi in [('png', 300), ('pdf', 300), ('svg', 300)]:
    fig.savefig(OUT / f'Figure4C_15genes.{ext}', dpi=dpi, facecolor='white')
fig.savefig(OUT / 'Figure4C_panel_preview.png', dpi=52.5, facecolor='white')
plt.close(fig)
validation = json.loads((ROOT / 'validation.json').read_text())
validation['figure'] = {'terms': len(rows), 'PNG_pixels': [3600, 3900 if tall else 3000], 'preview_pixels': [630, 682 if tall else 525], 'layout': args.layout, 'y_label_font_pt': 22 if tall else 18, 'bar_lengths': '-log10(Enrichr adjusted P-value)', 'n_label': 'Number of input genes overlapping each term', 'color_key': colors, 'plot_input': str(ROOT / 'Figure4C_plotted_terms.csv'), 'label_bounds_check': 'All axis labels, ticks and bar counts are within the canvas.'}
if tall:
    validation.pop('visual_check', None)
(OUT / 'validation.json').write_text(json.dumps(validation, indent=2) + '\n')
print('Saved PNG, PDF, SVG and 630-pixel preview.')
