# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Draw a readable Q1 enrichment panel including the endocrine and vascular axes."""
import os
import csv, json, os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/q1-matplotlib-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
OUT = Path(__file__).resolve().parent
rows = list(csv.DictReader((OUT / 'Figure2D_plotted_terms.csv').open()))
counts = {r['Keyword_category']: int(r['Significant_terms_n_ge_2']) for r in csv.DictReader((OUT / 'Q1_axis_counts.csv').open())}
palette = {'Neurodegenerative': '#598EAF', 'Inflammatory': '#B77798', 'Vascular': '#7B61A8', 'Endocrine': '#C59B55', 'Other': '#BAC3CA'}
ink = '#162E3E'
plt.rcParams.update({'font.family': 'Arial', 'text.color': ink, 'axes.labelcolor': ink, 'xtick.color': ink, 'ytick.color': ink, 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold', 'mathtext.default': 'regular'})
W, H = (13.5, 11.4)
fig = plt.figure(figsize=(W, H))
background = fig.add_axes([0, 0, 1, 1])
background.set(xlim=(0, W), ylim=(0, H))
background.axis('off')
heading = background.text(0.3, 10.98, 'Biological process enrichment', fontsize=26, fontweight='bold', va='center')
subheading = background.text(0.3, 10.48, '55 bridging genes · 16 representative GO terms', fontsize=21, va='center')
ax = fig.add_axes([0.565, 0.235, 0.4, 0.655])
cursor = 0.0
ys = []
for row in rows:
    height = 1.8 if '\n' in row['Plot_label'] else 1.0
    ys.append(cursor + height / 2)
    cursor += height
for row, y in zip(rows, ys):
    width = float(row['minus_log10_adjusted_P'])
    colors = row['Bar_colors'].split('; ')
    h = 0.62
    for i, color in enumerate(colors):
        ax.barh(y - h / 2 + h / len(colors) * (i + 0.5), width, height=h / len(colors), color=color, edgecolor='none')
    ax.add_patch(Rectangle((0, y - h / 2), width, h, facecolor='none', edgecolor='#71858C', lw=0.65))
    ax.text(width + 0.11, y, 'n=' + row['n_overlap'], fontsize=18, va='center')
ax.set_yticks(ys, [r['Plot_label'] for r in rows], fontsize=22)
ax.tick_params(axis='y', length=0, pad=12)
ax.tick_params(axis='x', labelsize=20, pad=7)
ax.set_ylim(cursor + 0.15, -0.15)
ax.set_xlim(0, 12)
ax.set_xticks([0, 3, 6, 9, 12])
ax.set_xlabel('$-\\log_{10}$(adjusted P value)', fontsize=23, fontweight='bold', labelpad=12)
ax.set_axisbelow(True)
ax.grid(axis='x', color='#E2E8EB', lw=0.7)
for s in ['top', 'right']:
    ax.spines[s].set_visible(False)
for s in ['left', 'bottom']:
    ax.spines[s].set_color('#A0AFB7')
    ax.spines[s].set_linewidth(0.8)
texts = [heading, subheading, *ax.get_yticklabels(), *ax.get_xticklabels(), ax.xaxis.label, *ax.texts]
texts.append(background.text(0.7, 1.69, 'Significant terms per keyword category', fontsize=20, va='center', fontweight='bold'))
for x, y, category in [(0.7, 1.14, 'Endocrine'), (4.2, 1.14, 'Inflammatory'), (8.35, 1.14, 'Neurodegenerative'), (0.7, 0.65, 'Vascular'), (4.2, 0.65, 'Other')]:
    background.add_patch(Rectangle((x, y - 0.11), 0.22, 0.22, facecolor=palette[category], edgecolor='#7E939F', lw=0.6))
    texts.append(background.text(x + 0.36, y, f'{category} ({counts[category]})', fontsize=20, va='center'))
texts.append(background.text(0.7, 0.21, 'n = overlapping genes; category counts include all 583 significant terms.', fontsize=17, va='center'))
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
for t in texts:
    box = t.get_window_extent(renderer)
    assert box.x0 >= 0 and box.y0 >= 0 and (box.x1 <= fig.bbox.width) and (box.y1 <= fig.bbox.height), (t.get_text(), box)
boxes = [t.get_window_extent(renderer) for t in ax.get_yticklabels()]
assert all((boxes[i].y0 > boxes[i + 1].y1 for i in range(len(boxes) - 1))), 'Overlapping term labels'
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(OUT / f'Figure2D_55genes_axes.{ext}', dpi=400, facecolor='white')
fig.savefig(OUT / 'Figure2D_preview.png', dpi=90, facecolor='white')
fig.savefig(OUT / 'Figure2D_panel_size.png', dpi=700 / W, facecolor='white')
plt.close(fig)
(OUT / 'plot_validation.json').write_text(json.dumps({'terms': len(rows), 'input_genes': 55, 'term_order': 'Enrichr adjusted P, raw P, term name', 'text_bounds': 'passed', 'y_label_overlap': 'none', 'bar_n_labels': 'Exact validated gene overlaps', 'bar_colors': 'Original nonexclusive keyword categories', 'multi_category_bar': 'Regulation of leukocyte adhesion to vascular endothelial cells: inflammatory + vascular, equal-height bands', 'legend_counts': counts, 'visual_review': 'Pending'}, indent=2) + '\n')
print('Saved Figure2D_55genes_axes PNG, PDF, SVG and previews.')
