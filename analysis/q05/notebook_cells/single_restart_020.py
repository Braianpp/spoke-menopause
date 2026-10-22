# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib_venn import venn3, venn3_circles
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica', 'Liberation Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['figure.dpi'] = 300
COLORS = {'apo': '#E63946', 'menopause': '#2A9D8F', 'cognitive': '#457B9D'}
OBSERVED_BAR_COLORS = {'apo_meno': '#D3C5C3', 'apo_cog': '#E9B6CA', 'meno_cog': '#9AE0E5', 'all_three': '#B1AEB8'}
try:
    df_null_summary
except NameError:
    null_summary_file = SAVE_PATH / f'degree_matched_null_summary{suffix}.csv'
    df_null_summary = pd.read_csv(null_summary_file)
    print(f'Loaded validation summary: {null_summary_file}')
try:
    jaccard_matrix
except NameError:
    jaccard_file = SAVE_PATH / f'jaccard_matrix_v4_compute{suffix}.csv'
    jaccard_matrix = pd.read_csv(jaccard_file, index_col=0)
    print(f'Loaded Jaccard matrix: {jaccard_file}')
validation_metrics = ['apo_meno', 'apo_cog', 'meno_cog', 'all_three']
metric_labels = {'apo_meno': 'APO–Meno', 'apo_cog': 'APO–Cog', 'meno_cog': 'Meno–Cog', 'all_three': 'All three'}
df_val = df_null_summary[df_null_summary['metric'].isin(validation_metrics)].copy()
df_val['label'] = df_val['metric'].map(metric_labels)
df_val['fold_over_null'] = df_val['observed'] / df_val['null_mean']
df_val['metric'] = pd.Categorical(df_val['metric'], categories=validation_metrics, ordered=True)
df_val = df_val.sort_values('metric').reset_index(drop=True)
jaccard_rows = [{'metric': 'apo_meno', 'label': 'APO–Meno', 'jaccard': float(jaccard_matrix.loc['APO', 'Menopause'])}, {'metric': 'apo_cog', 'label': 'APO–Cog', 'jaccard': float(jaccard_matrix.loc['APO', 'Cognitive'])}, {'metric': 'meno_cog', 'label': 'Meno–Cog', 'jaccard': float(jaccard_matrix.loc['Menopause', 'Cognitive'])}]
df_jaccard = pd.DataFrame(jaccard_rows)
fig_a, ax_a = plt.subplots(figsize=(6.2, 5.2))
venn = venn3([apo_genes, meno_genes, cog_genes], set_labels=('', '', ''), set_colors=(COLORS['apo'], COLORS['menopause'], COLORS['cognitive']), alpha=0.55, ax=ax_a)
venn3_circles([apo_genes, meno_genes, cog_genes], linewidth=2.5, color='white', ax=ax_a)
for lid in ['100', '010', '001', '110', '101', '011', '111']:
    label = venn.get_label_by_id(lid)
    if label:
        label.set_fontsize(13)
        label.set_fontweight('bold')
ax_a.set_title('A. Shared PPR-ranked genes', fontsize=14, fontweight='bold', loc='left', pad=20)
fig_a.tight_layout()
fig_a.savefig(f'${{SPOKE_DATA_ROOT}}/q7_figs/panel_A_shared_PPR_genes_venn{suffix}.svg', format='svg', bbox_inches='tight')
plt.show()
legend_handles_a = [Patch(facecolor=COLORS['apo'], edgecolor='none', alpha=0.55, label='APO'), Patch(facecolor=COLORS['menopause'], edgecolor='none', alpha=0.55, label='Menopause symptoms'), Patch(facecolor=COLORS['cognitive'], edgecolor='none', alpha=0.55, label='Cognitive / dementia')]
fig_leg_a, ax_leg_a = plt.subplots(figsize=(5.5, 0.6))
ax_leg_a.axis('off')
ax_leg_a.legend(handles=legend_handles_a, loc='center', ncol=3, frameon=False, fontsize=10)
fig_leg_a.savefig(f'${{SPOKE_DATA_ROOT}}/q7_figs/panel_A_legend{suffix}.svg', format='svg', bbox_inches='tight', transparent=True)
plt.show()
fig_b, ax_b = plt.subplots(figsize=(6.6, 5.2))
x = np.arange(len(df_val))
width = 0.35
ax_b.bar(x - width / 2, df_val['null_mean'], width, yerr=df_val['null_sd'], color='lightgray', edgecolor='gray', capsize=4)
ax_b.bar(x + width / 2, df_val['observed'], width, color=[OBSERVED_BAR_COLORS[m] for m in df_val['metric'].astype(str)], edgecolor='white', linewidth=1.5)
y_max = max(df_val['observed'].max(), df_val['null_mean'].max())
y_pad = y_max * 0.035
for i, (_, row) in enumerate(df_val.iterrows()):
    y = row['observed']
    p = row['empirical_p']
    fold = row['fold_over_null']
    p_text = 'p<0.001' if p <= 0.001 else f'p={p:.3f}'
    ax_b.text(i + width / 2, y + y_pad, f'{fold:.1f}×\n{p_text}', ha='center', va='bottom', fontsize=9, fontweight='bold')
ax_b.set_xticks(x)
ax_b.set_xticklabels(df_val['label'], fontsize=10)
ax_b.set_ylabel('Overlap genes', fontsize=11, fontweight='bold')
ax_b.set_title('B. Degree-matched permutation validation', fontsize=14, fontweight='bold', loc='left', pad=12)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.set_ylim(0, y_max * 1.3)
fig_b.tight_layout()
fig_b.savefig(f'${{SPOKE_DATA_ROOT}}/q7_figs/panel_B_degree_matched_validation{suffix}.svg', format='svg', bbox_inches='tight')
plt.show()
legend_handles_b = [Patch(facecolor='lightgray', edgecolor='gray', label='Degree-matched null'), Patch(facecolor=OBSERVED_BAR_COLORS['apo_meno'], edgecolor='none', label='Observed APO–Meno'), Patch(facecolor=OBSERVED_BAR_COLORS['apo_cog'], edgecolor='none', label='Observed APO–Cog'), Patch(facecolor=OBSERVED_BAR_COLORS['meno_cog'], edgecolor='none', label='Observed Meno–Cog'), Patch(facecolor=OBSERVED_BAR_COLORS['all_three'], edgecolor='none', label='Observed all three')]
fig_leg_b, ax_leg_b = plt.subplots(figsize=(7.0, 0.8))
ax_leg_b.axis('off')
ax_leg_b.legend(handles=legend_handles_b, loc='center', ncol=3, frameon=False, fontsize=9)
fig_leg_b.savefig(f'${{SPOKE_DATA_ROOT}}/q7_figs/panel_B_legend{suffix}.svg', format='svg', bbox_inches='tight', transparent=True)
plt.show()
fig_c, ax_c = plt.subplots(figsize=(5.6, 4.8))
x = np.arange(len(df_jaccard))
bars_c = ax_c.bar(x, df_jaccard['jaccard'], width=0.65, color=[OBSERVED_BAR_COLORS[m] for m in df_jaccard['metric']], edgecolor='white', linewidth=1.5)
for bar, value in zip(bars_c, df_jaccard['jaccard']):
    ax_c.text(bar.get_x() + bar.get_width() / 2, value + max(df_jaccard['jaccard']) * 0.035, f'{value:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
ax_c.set_xticks(x)
ax_c.set_xticklabels(df_jaccard['label'], fontsize=10)
ax_c.set_ylabel('Jaccard index', fontsize=11, fontweight='bold')
ax_c.set_title('C. Pairwise pathway similarity', fontsize=14, fontweight='bold', loc='left', pad=12)
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)
ax_c.set_ylim(0, max(df_jaccard['jaccard']) * 1.25)
fig_c.tight_layout()
fig_c.savefig(f'${{SPOKE_DATA_ROOT}}/q7_figs/panel_C_jaccard_similarity{suffix}.svg', format='svg', bbox_inches='tight')
plt.show()
legend_handles_c = [Patch(facecolor=OBSERVED_BAR_COLORS['apo_meno'], edgecolor='none', label='APO–Meno'), Patch(facecolor=OBSERVED_BAR_COLORS['apo_cog'], edgecolor='none', label='APO–Cog'), Patch(facecolor=OBSERVED_BAR_COLORS['meno_cog'], edgecolor='none', label='Meno–Cog')]
fig_leg_c, ax_leg_c = plt.subplots(figsize=(4.8, 0.6))
ax_leg_c.axis('off')
ax_leg_c.legend(handles=legend_handles_c, loc='center', ncol=3, frameon=False, fontsize=9)
fig_leg_c.savefig(f'${{SPOKE_DATA_ROOT}}/q7_figs/panel_C_legend{suffix}.svg', format='svg', bbox_inches='tight', transparent=True)
plt.show()
print('Saved SVG files:')
print(f'- panel_A_shared_PPR_genes_venn{suffix}.svg')
print(f'- panel_A_legend{suffix}.svg')
print(f'- panel_B_degree_matched_validation{suffix}.svg')
print(f'- panel_B_legend{suffix}.svg')
print(f'- panel_C_jaccard_similarity{suffix}.svg')
print(f'- panel_C_legend{suffix}.svg')
