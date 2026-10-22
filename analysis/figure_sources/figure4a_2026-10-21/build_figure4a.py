# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Correct Figure 4A highlight selection using the seven curated MeSH IDs.

Uses saved dementia-seeded PageRank results. No graph query, PageRank rerun,
network change, or enrichment rerun is performed.
"""
import os
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import os
OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--legend-inside', action='store_true', help='Write an additional PNG with the legend inside the plot.')
args = parser.parse_args()
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/q3-matplotlib-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.stats import gaussian_kde
SOURCE = ROOT / 'audit_2026-09-30/excel_stage/Q3_SymptomPPR_Dementia.csv'
CURATED = ROOT / 'audit_2026-09-30/excel_stage/Q3_CommonSymptoms_Summary.csv'
PAYLOAD = ROOT / 'supplementary_tables_2026-10-01/payload_submission.json'

def read_csv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))
rows = read_csv(SOURCE)
by_id = {r['mesh_id']: r for r in rows}
curated = read_csv(CURATED)
selected_ids = ['D010292', 'D059352', 'D007319', 'D020069', 'D050177', 'D019584', 'D007039']
assert len(rows) == len(by_id) == 1891
assert set(selected_ids) == {r['mesh_id'] for r in curated}
assert set(selected_ids) == {r['mesh_id'] for r in rows if r['common_perimenopausal_symptom (7)'] == 'True'}
assert 'D000068376' not in selected_ids
assert by_id['D000068376']['symptom'] == 'Compassion Fatigue'
payload = json.loads(PAYLOAD.read_text())
table = next((t for t in payload['tables'] if t['number'] == 13))
selected_block = next((b for b in table['blocks'] if b['title'] == 'Selected perimenopausal symptoms'))
background_block = next((b for b in table['blocks'] if b['title'] == 'All symptom nodes'))
assert {r[1] for r in selected_block['rows']} == set(selected_ids)
assert len(background_block['rows']) == 1891
for r in background_block['rows']:
    assert math.isclose(float(r[2]), float(by_id[r[1]]['pagerank']), rel_tol=1e-12, abs_tol=0)
scores = np.array([float(r['pagerank']) for r in rows])
transformed = np.log10(scores + 1e-10)
visible = transformed > -9
values = transformed[visible]
assert visible.sum() == 1398
assert int((scores == 0).sum()) == 493
assert np.array_equal(visible, scores > 0)
assert visible[rows.index(by_id['D000068376'])]
counts, edges = np.histogram(values, bins=50)
assert counts.sum() == 1398 and counts.max() == 56
assert counts[:13].tolist() == [1, 1, 5, 0, 1, 6, 45, 39, 4, 26, 45, 30, 40]
ink = '#162E3E'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 17, 'text.color': ink, 'axes.labelcolor': ink, 'xtick.color': ink, 'ytick.color': ink, 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
colors = {'D010292': '#0072B2', 'D059352': '#D55E00', 'D007319': '#9467BD', 'D020069': '#009E73', 'D050177': '#B19A2E', 'D019584': '#CC79A7', 'D007039': '#7E8B94'}
names = {'D059352': 'Musculoskeletal pain', 'D020069': 'Shoulder pain', 'D019584': 'Hot flashes'}
fig, ax = plt.subplots(figsize=(12, 8.4))
fig.subplots_adjust(left=0.105, right=0.978, bottom=0.145 if args.legend_inside else 0.32, top=0.97)
ax.hist(values, bins=edges, color='#D9D9D9', edgecolor='white', linewidth=0.65, zorder=1)
grid = np.linspace(values.min(), values.max(), 200)
ax.plot(grid, gaussian_kde(values)(grid) * len(values) * (edges[1] - edges[0]), color='#B8BFC3', lw=2.2, zorder=2)
handles = []
plotted = []
for mid in selected_ids:
    r = by_id[mid]
    score = float(r['pagerank'])
    name = names.get(mid, r['symptom'])
    x = math.log10(score + 1e-10)
    if score > 0:
        ax.axvline(x, color=colors[mid], linestyle='--', lw=2.8, zorder=4)
        handles.append(Line2D([], [], color=colors[mid], ls='--', lw=2.8, label=name))
        shown = True
    else:
        shown = False
    plotted.append({'symptom': name, 'mesh_id': mid, 'pagerank': score, 'rank_of_1891_symptoms': int(r['rank_of_all_symptoms']), 'log10_pagerank_plus_1e-10': x, 'vertical_line_shown': shown, 'color': colors[mid]})
assert sum((r['vertical_line_shown'] for r in plotted)) == 6
ax.set_xlim(-9, -4)
ax.set_ylim(0, 60)
ax.set_xticks([-9, -8, -7, -6, -5, -4])
ax.set_yticks([0, 10, 20, 30, 40, 50, 60])
ax.set_xlabel('$\\log_{10}$(PageRank + $10^{-10}$)', fontsize=19, fontweight='bold', labelpad=12)
ax.set_ylabel('Symptom count', fontsize=19, fontweight='bold', labelpad=12)
ax.tick_params(labelsize=17, length=4, pad=5)
ax.set_axisbelow(True)
ax.grid(color='#DCE3E6', linewidth=0.7)
for spine in ['left', 'bottom']:
    ax.spines[spine].set_color('#9EAEB4')
    ax.spines[spine].set_linewidth(0.8)
if args.legend_inside:
    legend = ax.legend(handles=handles, title='Selected symptoms', loc='upper left', bbox_to_anchor=(0.012, 0.993), frameon=True, facecolor='white', edgecolor='none', framealpha=0.9, fontsize=17, title_fontsize=18, handlelength=2.0, borderpad=0.5, labelspacing=0.5)
else:
    legend = fig.legend(handles=handles, title='Selected perimenopausal symptoms', loc='lower center', bbox_to_anchor=(0.5415, 0.017), frameon=False, ncol=2, fontsize=17, title_fontsize=18, handlelength=2.0, borderpad=0, labelspacing=0.5, columnspacing=1.6)
legend.get_title().set_fontweight('bold')
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
for artist in [ax.xaxis.label, ax.yaxis.label, *ax.get_xticklabels(), *ax.get_yticklabels(), legend]:
    box = artist.get_window_extent(renderer)
    assert box.x0 >= 0 and box.y0 >= 0 and (box.x1 <= fig.bbox.width) and (box.y1 <= fig.bbox.height)
assert 'Compassion Fatigue' not in [t.get_text() for t in legend.get_texts()]
assert all(('Sparse hair' not in t.get_text() for t in legend.get_texts()))
assert len(legend.get_texts()) == 6
stem = 'Figure4A_legend_inside' if args.legend_inside else 'Figure4A_selected_symptoms'
for extension in ['png'] if args.legend_inside else ['png', 'pdf', 'svg']:
    fig.savefig(OUT / f'{stem}.{extension}', dpi=600, facecolor='white')
preview = 'Figure4A_legend_inside_preview.png' if args.legend_inside else 'Figure4A_preview.png'
fig.savefig(OUT / preview, dpi=90, facecolor='white')
plt.close(fig)
if args.legend_inside:
    check = {'figure': stem + '.png', 'format': 'PNG', 'dpi': 600, 'pixels': [7200, 5040], 'legend_placement': 'Inside axes, upper left', 'legend_entries': 6, 'displayed_ids': [r['mesh_id'] for r in plotted if r['vertical_line_shown']], 'background_symptoms': int(visible.sum()), 'histogram_counts': counts.tolist(), 'label_bounds_check': 'passed', 'visual_review': 'pending'}
    (OUT / 'validation_legend_inside.json').write_text(json.dumps(check, indent=2) + '\n')
    print('Saved additional PNG with legend inside; existing exports preserved.')
    raise SystemExit(0)
with (OUT / 'Figure4A_selected_symptoms.csv').open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(plotted[0]))
    writer.writeheader()
    writer.writerows(plotted)
with (OUT / 'Figure4A_histogram_bins.csv').open('w', newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(['bin_lower', 'bin_upper', 'symptom_count'])
    writer.writerows(zip(edges[:-1], edges[1:], counts))
caption = 'a, Distribution of dementia-seeded personalized PageRank across the 1,398 SPOKE symptom nodes with nonzero scores, shown as log10(PageRank + 10^-10). Colored dashed lines mark the six selected perimenopausal symptoms with nonzero scores: paresthesia, musculoskeletal pain, insomnia, shoulder pain, overweight and hot flashes. Selected symptoms with zero scores are omitted from the display. Of all 1,891 symptom nodes, 493 have zero scores and are outside this range. The gray curve is a kernel density estimate scaled to the histogram counts. Selected symptoms were matched by their exact MeSH identifiers. Scores, ranks and identifiers are provided in Supplementary Table 13.\n'
(OUT / 'Figure4A_caption.txt').write_text(caption)
(OUT / 'Figure4A_manuscript_update.txt').write_text('Replace Figure 4A and its panel-a caption with the corrected panel and accompanying caption.\n\nMethods clarification: Symptoms highlighted in Figure 4a were selected by exact MeSH identifier from the curated seven-symptom list (Supplementary Table 13); the figure displays the six selected symptoms with nonzero PageRank scores.\n\nThe Q3 Results scores and ranks for paresthesia and insomnia do not change. The 15-gene Figure 4B network and Figure 4C enrichment do not change. Supplementary Table 13 already flags Compassion Fatigue as unselected; it remains in the all-symptom background. No workbook change is required.\nSparse hair has a zero score and is omitted from the figure and legend at Braian’s request. It is retained in the complete supplementary results.\n')
validation = {'selection_rule': 'Exact membership in the seven curated MeSH IDs; display only PageRank > 0', 'selected_symptoms': plotted, 'all_symptom_nodes': len(rows), 'background_with_nonzero_score': int(visible.sum()), 'zero_scores': int((scores == 0).sum()), 'histogram_bins': 50, 'histogram_maximum_count': int(counts.max()), 'histogram_matches_original_bins': True, 'submission_table_13_scores_and_selection_match': True, 'compassion_fatigue_highlighted': False, 'compassion_fatigue_retained_in_background': True, 'legend_entries': 6, 'zero_score_selected_symptoms_omitted_from_figure': True, 'pagerank_rerun': False, 'network_changed': False, 'enrichment_changed': False, 'font_family': 'Arial', 'axis_label_font_pt': 19, 'axis_label_font_weight': 'bold', 'text_color': ink, 'png_dpi': 600, 'png_dimensions': [7200, 5040], 'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [SOURCE, CURATED, PAYLOAD]}, 'label_bounds_check': 'passed', 'legend_placement': 'Below the axes; no histogram bars obscured', 'visual_review': 'pending'}
(OUT / 'validation.json').write_text(json.dumps(validation, indent=2) + '\n')
print(json.dumps({k: v for k, v in validation.items() if k not in ['source_sha256', 'selected_symptoms']}, indent=2))
