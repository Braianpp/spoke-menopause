# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""
DWPC Visualizations for MHT → Dementia
=======================================
Run after the main MHT analysis. Reads dwpc CSV and creates presentation-ready plots.
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings('ignore')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica', 'Liberation Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['figure.dpi'] = 300
SAVE_PATH = '${SPOKE_DATA_ROOT}/'
COLORS = {'Conjugated estrogens': '#E63946', 'Estradiol': '#F4A261', 'Progesterone': '#2A9D8F', 'Medroxyprogesterone acetate': '#264653', 'Norethindrone acetate': '#E76F51'}
df = pd.read_csv(f'{SAVE_PATH}mht_dementia_dwpc_2hops.csv')
print(f'Loaded {len(df):,} DWPC results')
print(f"Compounds: {df['source_name'].unique()}")
print(f"Targets: {df['target_name'].unique()}")
df_alz = df[df['target_name'] == "Alzheimer's disease"].copy()
top_per_compound = []
for compound in df_alz['source_name'].unique():
    sub = df_alz[df_alz['source_name'] == compound].nlargest(3, 'dwpc')
    top_per_compound.append(sub)
df_top = pd.concat(top_per_compound).sort_values('dwpc', ascending=True)

def clean_metapath(mp):
    mp = mp.replace('TREATS_CtD', 'Treats')
    mp = mp.replace('ASSOCIATES_VaP', 'VariantAssoc')
    mp = mp.replace('MEASURES_CLmC', 'LabMeasure')
    mp = mp.replace('UPREGULATES_CuG', 'Upreg')
    mp = mp.replace('DOWNREGULATES_CdG', 'Downreg')
    mp = mp.replace('INTERACTS_PiC', 'ProteinInt')
    mp = mp.replace('ISA_DiD', 'IsA')
    mp = mp.replace('INCREASEDIN_PiD', 'IncreasedIn')
    mp = mp.replace('PARTOF_CpPG', 'PartOf')
    mp = mp.replace(' → ', '→')
    return mp
df_top['metapath_clean'] = df_top['metapath'].apply(clean_metapath)
df_top['label'] = df_top['source_name'] + '\n' + df_top['metapath_clean']
fig1, ax1 = plt.subplots(figsize=(14, 10))
bar_colors = [COLORS.get(c, '#888888') for c in df_top['source_name']]
bars = ax1.barh(range(len(df_top)), df_top['dwpc'], color=bar_colors, edgecolor='white', linewidth=0.8, alpha=0.85)
ax1.set_yticks(range(len(df_top)))
ax1.set_yticklabels(df_top['label'], fontsize=9)
ax1.set_xlabel('DWPC Score', fontsize=12, fontweight='bold')
ax1.set_title("Top Metapaths: MHT Compounds → Alzheimer's Disease\n(Top 3 per compound, DWPC-ranked)", fontsize=14, fontweight='bold')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
for bar, (_, row) in zip(bars, df_top.iterrows()):
    w = bar.get_width()
    ax1.text(w + 0.001, bar.get_y() + bar.get_height() / 2, f"{row['path_count']:,} paths", va='center', fontsize=8, color='gray')
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=COLORS[c], label=c) for c in COLORS]
ax1.legend(handles=legend_elements, loc='lower right', fontsize=9, framealpha=0.9)
plt.tight_layout()
fig1.savefig(f'{SAVE_PATH}dwpc_top_metapaths_alzheimers.png', dpi=300, bbox_inches='tight')
fig1.savefig(f'{SAVE_PATH}dwpc_top_metapaths_alzheimers.pdf', dpi=300, bbox_inches='tight')
print('✅ Figure 1: Top metapaths bar chart saved')
plt.show()
best = df.groupby(['source_name', 'target_name']).agg(best_dwpc=('dwpc', 'max'), total_paths=('path_count', 'sum'), n_active_metapaths=('dwpc', lambda x: (x > 0).sum())).reset_index()
fig2, ax2 = plt.subplots(figsize=(12, 7))
compound_order = best.groupby('source_name')['best_dwpc'].max().sort_values(ascending=False).index.tolist()
target_order = ["Alzheimer's disease", 'Vascular dementia', 'Dementia']
for i, compound in enumerate(compound_order):
    sub = best[best['source_name'] == compound]
    for _, row in sub.iterrows():
        j = target_order.index(row['target_name'])
        size = max(row['total_paths'] / 50, 30)
        ax2.scatter(j, i, s=min(size, 800), c=COLORS.get(compound, '#888'), alpha=0.7, edgecolors='white', linewidth=1.5)
        ax2.text(j, i + 0.25, f"{row['best_dwpc']:.3f}", ha='center', va='bottom', fontsize=9, fontweight='bold')
ax2.set_yticks(range(len(compound_order)))
ax2.set_yticklabels(compound_order, fontsize=11)
ax2.set_xticks(range(len(target_order)))
ax2.set_xticklabels(target_order, fontsize=11)
ax2.set_title('Best DWPC Score: MHT Compound × Dementia Target\n(Bubble size ∝ total path count)', fontsize=14, fontweight='bold')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_xlim(-0.5, 2.5)
ax2.set_ylim(-0.5, len(compound_order) - 0.5)
ax2.invert_yaxis()
plt.tight_layout()
fig2.savefig(f'{SAVE_PATH}dwpc_bubble_compound_target.png', dpi=300, bbox_inches='tight')
fig2.savefig(f'{SAVE_PATH}dwpc_bubble_compound_target.pdf', dpi=300, bbox_inches='tight')
print('✅ Figure 2: Bubble chart saved')
plt.show()

def classify_metapath(mp):
    if 'ASSOCIATES_VaP' in mp:
        return 'Genetic variant\nassociation'
    elif 'UPREGULATES' in mp or 'DOWNREGULATES' in mp:
        return 'Gene\nregulation'
    elif 'TREATS_CtD → TREATS_CtD' in mp or (mp.count('TREATS_CtD') >= 2 and 'ISA' not in mp):
        return 'Therapeutic\noverlap'
    elif 'INTERACTS' in mp:
        return 'Protein\ninteraction'
    elif 'MEASURES' in mp:
        return 'Biomarker\noverlap'
    elif 'ISA' in mp:
        return 'Disease\nhierarchy'
    else:
        return 'Other'
df_alz['mechanism'] = df_alz['metapath'].apply(classify_metapath)
mech_summary = df_alz.groupby(['source_name', 'mechanism'])['dwpc'].sum().reset_index()
mech_pivot = mech_summary.pivot(index='mechanism', columns='source_name', values='dwpc').fillna(0)
mech_order = mech_pivot.sum(axis=1).sort_values(ascending=False).index.tolist()
mech_pivot = mech_pivot.loc[mech_order]
compound_order_fig3 = ['Medroxyprogesterone acetate', 'Progesterone', 'Estradiol', 'Conjugated estrogens', 'Norethindrone acetate']
mech_pivot = mech_pivot[[c for c in compound_order_fig3 if c in mech_pivot.columns]]
fig3, ax3 = plt.subplots(figsize=(14, 7))
x = np.arange(len(mech_pivot))
width = 0.15
n_compounds = len(mech_pivot.columns)
for i, compound in enumerate(mech_pivot.columns):
    offset = (i - n_compounds / 2 + 0.5) * width
    bars = ax3.bar(x + offset, mech_pivot[compound], width, label=compound, color=COLORS.get(compound, '#888'), edgecolor='white', linewidth=0.5, alpha=0.85)
ax3.set_xticks(x)
ax3.set_xticklabels(mech_pivot.index, fontsize=10)
ax3.set_ylabel('Cumulative DWPC Score', fontsize=12, fontweight='bold')
ax3.set_title("MHT → Alzheimer's: Connectivity by Mechanism Type\n(Sum of DWPC scores across all metapaths of each type)", fontsize=14, fontweight='bold')
ax3.legend(fontsize=8, loc='upper right', framealpha=0.9)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
plt.tight_layout()
fig3.savefig(f'{SAVE_PATH}dwpc_mechanism_comparison.png', dpi=300, bbox_inches='tight')
fig3.savefig(f'{SAVE_PATH}dwpc_mechanism_comparison.pdf', dpi=300, bbox_inches='tight')
print('✅ Figure 3: Mechanism comparison saved')
plt.show()
fig4 = plt.figure(figsize=(22, 8))
gs = GridSpec(1, 3, figure=fig4, wspace=0.35)
ax4a = fig4.add_subplot(gs[0, 0])
best_pivot = best.pivot(index='source_name', columns='target_name', values='best_dwpc')
best_pivot = best_pivot[target_order]
row_order = best_pivot["Alzheimer's disease"].sort_values(ascending=False).index
best_pivot = best_pivot.loc[row_order]
sns.heatmap(best_pivot, annot=True, fmt='.3f', cmap='Purples', ax=ax4a, linewidths=0.5, linecolor='white', cbar_kws={'label': 'DWPC', 'shrink': 0.8})
ax4a.set_title('A. Best DWPC per\nCompound–Target Pair', fontsize=13, fontweight='bold')
ax4a.set_ylabel('')
ax4a.set_xlabel('')
ax4b = fig4.add_subplot(gs[0, 1])
top2 = []
for compound in compound_order_fig3:
    sub = df_alz[df_alz['source_name'] == compound].nlargest(2, 'dwpc')
    top2.append(sub)
df_top2 = pd.concat(top2).sort_values('dwpc', ascending=True)
df_top2['metapath_clean'] = df_top2['metapath'].apply(clean_metapath)
df_top2['label'] = df_top2['source_name'].apply(lambda x: x[:15]) + '\n' + df_top2['metapath_clean']
bar_colors2 = [COLORS.get(c, '#888') for c in df_top2['source_name']]
ax4b.barh(range(len(df_top2)), df_top2['dwpc'], color=bar_colors2, edgecolor='white', linewidth=0.8, alpha=0.85)
ax4b.set_yticks(range(len(df_top2)))
ax4b.set_yticklabels(df_top2['label'], fontsize=8)
ax4b.set_xlabel('DWPC Score', fontsize=10)
ax4b.set_title("B. Top Metapaths → Alzheimer's\n(Top 2 per compound)", fontsize=13, fontweight='bold')
ax4b.spines['top'].set_visible(False)
ax4b.spines['right'].set_visible(False)
ax4c = fig4.add_subplot(gs[0, 2])
mech_for_stack = df_alz.groupby(['source_name', 'mechanism'])['dwpc'].sum().reset_index()
mech_for_stack_pivot = mech_for_stack.pivot(index='source_name', columns='mechanism', values='dwpc').fillna(0)
mech_for_stack_pivot['total'] = mech_for_stack_pivot.sum(axis=1)
mech_for_stack_pivot = mech_for_stack_pivot.sort_values('total', ascending=True)
mech_for_stack_pivot = mech_for_stack_pivot.drop('total', axis=1)
mech_colors = {'Genetic variant\nassociation': '#8E44AD', 'Gene\nregulation': '#2ECC71', 'Therapeutic\noverlap': '#3498DB', 'Protein\ninteraction': '#E67E22', 'Biomarker\noverlap': '#1ABC9C', 'Disease\nhierarchy': '#95A5A6', 'Other': '#BDC3C7'}
col_order = mech_for_stack_pivot.sum().sort_values(ascending=False).index
mech_for_stack_pivot = mech_for_stack_pivot[col_order]
left = np.zeros(len(mech_for_stack_pivot))
for mech in mech_for_stack_pivot.columns:
    vals = mech_for_stack_pivot[mech].values
    mech_label = mech.replace('\n', ' ')
    ax4c.barh(range(len(mech_for_stack_pivot)), vals, left=left, color=mech_colors.get(mech, '#BDC3C7'), label=mech_label, edgecolor='white', linewidth=0.5)
    left += vals
ax4c.set_yticks(range(len(mech_for_stack_pivot)))
ax4c.set_yticklabels(mech_for_stack_pivot.index, fontsize=15)
ax4c.set_xlabel('Total DWPC (sum across metapaths)', fontsize=10)
ax4c.set_title("C. Mechanism Breakdown\n(Alzheimer's)", fontsize=13, fontweight='bold')
ax4c.legend(fontsize=7, loc='lower right', framealpha=0.9)
ax4c.spines['top'].set_visible(False)
ax4c.spines['right'].set_visible(False)
fig4.suptitle('MHT → Dementia: DWPC Metapath Analysis', fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
fig4.savefig(f'{SAVE_PATH}dwpc_presentation_panel.png', dpi=300, bbox_inches='tight')
fig4.savefig(f'{SAVE_PATH}dwpc_presentation_panel.pdf', dpi=300, bbox_inches='tight')
print('✅ Figure 4: Combined presentation panel saved')
plt.show()
print('\n' + '=' * 80)
print('FIGURES SAVED')
print('=' * 80)
print(f'  1. dwpc_top_metapaths_alzheimers.png/pdf  — Top 3 metapaths per compound')
print(f'  2. dwpc_bubble_compound_target.png/pdf     — Bubble: compound × target')
print(f'  3. dwpc_mechanism_comparison.png/pdf       — Grouped bars by mechanism type')
print(f'  4. dwpc_presentation_panel.png/pdf         — Combined 3-panel (best for slides)')
