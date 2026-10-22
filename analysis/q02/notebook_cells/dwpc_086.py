# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import rcParams
from neo4j import GraphDatabase, basic_auth
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
rcParams['axes.linewidth'] = 0.6
SAVE_PREFIX = 'pm_all_diseases'
SELECTED_TARGETS = {'DOID:8725': 'Vascular dementia', 'DOID:9255': 'Frontotemporal dementia', 'DOID:14330': "Parkinson's disease"}
OTHER_NEURO = {'DOID:10652': "Alzheimer's disease", 'DOID:332': 'Amyotrophic lateral sclerosis', 'DOID:12858': "Huntington's disease", 'DOID:678': 'Multiple sclerosis', 'DOID:1289': 'Neurodegenerative disease', 'DOID:12217': 'Lewy body dementia', 'DOID:2474': 'Dementia'}
C_SELECTED = {'Vascular dementia': '#E74C3C', 'Frontotemporal dementia': '#3498DB', "Parkinson's disease": '#2ECC71'}
C_OTHER_NEURO = '#9B59B6'
C_HIST = '#D5D5D5'
C_TEXT = '#333333'
C_SUB = '#777777'
print('=' * 70)
print('STEP 1: Querying 3-hop network proximity from PM to all diseases')
print('=' * 70)
query = "\nMATCH (pm:Disease {identifier: 'DOID:10787'})-[:ASSOCIATES_DaG]-(g1:Gene)-[]-(g2:Gene)-[:ASSOCIATES_DaG]-(d:Disease)\nWHERE d <> pm AND g1 <> g2\nWITH d.identifier AS disease_id, d.name AS disease,\n     count(DISTINCT g2) AS shared_3hop\nRETURN disease_id, disease, shared_3hop\nORDER BY shared_3hop DESC\n"
print('  Running query (this takes ~30-60 seconds)...')
with open_spoke(SPOKE_User, SPOKE_Password, SPOKE_DEV_URI) as db:
    session = db.session()
    dwpc_all_df = pd.DataFrame(session.run(query).data())
print(f'  ✅ Done: {len(dwpc_all_df)} diseases retrieved')
dwpc_all_df = dwpc_all_df.sort_values('shared_3hop', ascending=False).reset_index(drop=True)
dwpc_all_df['rank'] = range(1, len(dwpc_all_df) + 1)
dwpc_all_df.to_csv(f'{SAVE_PREFIX}_3hop.csv', index=False)
print(f'  Saved: {SAVE_PREFIX}_3hop.csv')
n_total = len(dwpc_all_df)
n_nonzero = (dwpc_all_df['shared_3hop'] > 0).sum()
print(f'  Diseases with connectivity > 0: {n_nonzero} / {n_total}')
print(f"\n{'=' * 70}")
print('STEP 2: Disease rankings')
print('=' * 70)
print(f'\n  Selected targets:')
for did, dname in SELECTED_TARGETS.items():
    row = dwpc_all_df[dwpc_all_df['disease_id'] == did]
    if len(row) > 0:
        r = row.iloc[0]
        pct = r['rank'] / n_total * 100
        print(f"    {dname:<30s}  Shared={int(r['shared_3hop'])}  Rank {int(r['rank'])}/{n_total} (top {pct:.1f}%)")
    else:
        print(f'    {dname:<30s}  NOT FOUND')
print(f'\n  Other neurodegenerative diseases:')
for did, dname in OTHER_NEURO.items():
    row = dwpc_all_df[dwpc_all_df['disease_id'] == did]
    if len(row) > 0:
        r = row.iloc[0]
        pct = r['rank'] / n_total * 100
        print(f"    {dname:<30s}  Shared={int(r['shared_3hop'])}  Rank {int(r['rank'])}/{n_total} (top {pct:.1f}%)")
    else:
        print(f'    {dname:<30s}  NOT FOUND')
print(f'\n  Top 20 diseases overall:')
for _, r in dwpc_all_df.head(20).iterrows():
    flag = ''
    if r['disease_id'] in SELECTED_TARGETS:
        flag = ' ★'
    elif r['disease_id'] in OTHER_NEURO:
        flag = ' ●'
    print(f"    #{int(r['rank']):>4d}  {r['disease']:<45s}  n={int(r['shared_3hop'])}{flag}")
print(f"\n{'=' * 70}")
print('STEP 3: Generating distribution figure')
print('=' * 70)
positive = dwpc_all_df[dwpc_all_df['shared_3hop'] > 0].copy()
fig, ax = plt.subplots(figsize=(11, 5.5))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')
ax.set_axisbelow(True)
ax.yaxis.grid(True, color='#F0F0F0', linewidth=0.4)
ax.hist(positive['shared_3hop'], bins=60, color=C_HIST, edgecolor='white', linewidth=0.3, alpha=0.85, zorder=2)
y_max = ax.get_ylim()[1]
next_y = y_max * 0.95
for did, dname in SELECTED_TARGETS.items():
    row = dwpc_all_df[dwpc_all_df['disease_id'] == did]
    if len(row) == 0:
        continue
    score = row.iloc[0]['shared_3hop']
    rank = int(row.iloc[0]['rank'])
    color = C_SELECTED[dname]
    ax.axvline(score, color=color, linewidth=2.5, linestyle='--', zorder=4, alpha=0.9)
    ax.text(score, next_y, f'  {dname}\n  (#{rank}, n={int(score)})', fontsize=7.5, fontweight='bold', color=color, va='top', ha='left', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=color, alpha=0.9, linewidth=0.8), zorder=5)
    next_y -= y_max * 0.15
neuro_found = []
for did, dname in OTHER_NEURO.items():
    row = dwpc_all_df[dwpc_all_df['disease_id'] == did]
    if len(row) == 0:
        continue
    score = row.iloc[0]['shared_3hop']
    if score > 0:
        ax.axvline(score, color=C_OTHER_NEURO, linewidth=1.2, linestyle='-', zorder=3, alpha=0.5)
        neuro_found.append((dname, score, int(row.iloc[0]['rank'])))
neuro_found.sort(key=lambda x: x[1], reverse=True)
label_y = y_max * 0.42
for dname, score, rank in neuro_found[:5]:
    ax.text(score, label_y, f'  {dname} (#{rank})', fontsize=6.5, color=C_OTHER_NEURO, va='top', ha='left', fontstyle='italic', alpha=0.8, zorder=5)
    label_y -= y_max * 0.06
ax.set_xlim(0, 200)
ax.set_xlabel('Shared 3-Hop Gene Intermediates (PM → Gene → Gene → Disease)', fontsize=10, fontweight='bold', labelpad=6)
ax.set_ylabel('Number of Diseases', fontsize=10, fontweight='bold', labelpad=6)
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)
ax.set_title('Network Proximity: Premature Menopause to All Diseases in SPOKE', fontsize=12, fontweight='bold', loc='left', pad=10)
ax.text(0, 1.02, f'', transform=ax.transAxes, fontsize=7.5, color=C_SUB, va='bottom')
legend_elements = [Patch(facecolor=C_HIST, edgecolor='white', label='All diseases')]
for dname in SELECTED_TARGETS.values():
    legend_elements.append(Line2D([0], [0], color=C_SELECTED[dname], linewidth=2.5, linestyle='--', label=dname))
legend_elements.append(Line2D([0], [0], color=C_OTHER_NEURO, linewidth=1.2, linestyle='-', alpha=0.6, label='Other neurodegenerative'))
ax.legend(handles=legend_elements, loc='upper right', fontsize=7.5, frameon=True, framealpha=0.95, edgecolor='#DDD', fancybox=True)
stats_text = f"Diseases: {n_total:,}\nConnected: {n_nonzero}\nMedian: {positive['shared_3hop'].median():.0f}\nMean: {positive['shared_3hop'].mean():.1f}"
ax.text(0.98, 0.72, stats_text, transform=ax.transAxes, fontsize=7, va='top', ha='right', color=C_TEXT, bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#DDD', alpha=0.95))
plt.tight_layout()
fig.savefig(f'{SAVE_PREFIX}_3hop_distribution.png', dpi=300, bbox_inches='tight', facecolor='white')
print(f'  ✅ Saved: {SAVE_PREFIX}_3hop_distribution.png')
plt.show()
print(f"\n{'=' * 70}")
print('Top 30 diseases by 3-hop proximity to Premature Menopause')
print('=' * 70)
top30 = dwpc_all_df.head(30).copy()
top30['category'] = 'Other'
for did in SELECTED_TARGETS:
    top30.loc[top30['disease_id'] == did, 'category'] = '★ Selected'
for did in OTHER_NEURO:
    top30.loc[top30['disease_id'] == did, 'category'] = '● Neurodegen.'
print(top30[['rank', 'disease', 'shared_3hop', 'category']].to_string(index=False))
print(f"\n{'=' * 70}")
print('SUMMARY')
print('=' * 70)
print(f'\nMetric: 3-hop shared gene intermediates\nPath:   PM -[ASSOCIATES_DaG]- Gene1 -[any]- Gene2 -[ASSOCIATES_DaG]- Disease\nThis matches the length-3 metapaths used in the DWPC analysis.\n\nTotal diseases assessed:  {n_total:,}\nWith connectivity > 0:    {n_nonzero}\n\nSelected targets:')
for did, dname in SELECTED_TARGETS.items():
    row = dwpc_all_df[dwpc_all_df['disease_id'] == did]
    if len(row) > 0:
        r = row.iloc[0]
        pct = r['rank'] / n_total * 100
        print(f"  {dname:<30s}  n={int(r['shared_3hop'])}  #{int(r['rank'])} / {n_total} (top {pct:.1f}%)")
print(f'\nOther neurodegenerative diseases:')
for did, dname in OTHER_NEURO.items():
    row = dwpc_all_df[dwpc_all_df['disease_id'] == did]
    if len(row) > 0:
        r = row.iloc[0]
        pct = r['rank'] / n_total * 100
        print(f"  {dname:<30s}  n={int(r['shared_3hop'])}  #{int(r['rank'])} / {n_total} (top {pct:.1f}%)")
print(f'\nFiles saved:\n  • {SAVE_PREFIX}_3hop.csv              — Full proximity scores for all diseases\n  • {SAVE_PREFIX}_3hop_distribution.png — Distribution histogram\n')
