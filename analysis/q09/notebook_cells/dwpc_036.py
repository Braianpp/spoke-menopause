# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""
MHT → Dementia: Compound-Disease Association Analysis
======================================================
Research Question: Are menopausal hormone therapy (MHT) compounds more
associated with dementia than random compounds in SPOKE?

Approach:
  1. Measure MHT–dementia proximity via BFS gene overlap + PPR embedding
  2. Compute DWPC (degree-weighted path counts) for MHT→Dementia metapaths
  3. Build null distribution from random Compound nodes
  4. Per-compound breakdown: which MHT is most connected to dementia?

Requires: G (SPOKE networkx graph), NetworkDWPC calculator
"""
import os
from scipy.stats import fisher_exact
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
from scipy import stats
import random
import time
import warnings
warnings.filterwarnings('ignore')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica', 'Liberation Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['figure.dpi'] = 300
SAVE_PATH = '${SPOKE_DATA_ROOT}/'
MAX_HOPS = 2
N_PERMUTATIONS = 1000
RANDOM_SEED = 42
USE_DWPC = True
DWPC_MAX_LENGTH = 3
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)
COLORS = {'mht': '#8E44AD', 'random': '#BDC3C7', 'dementia': '#457B9D', 'significant': '#E74C3C', 'conj_estrogens': '#E63946', 'estradiol': '#F4A261', 'progesterone': '#2A9D8F', 'medroxy': '#264653', 'norethindrone': '#E76F51'}
mht_compounds = {'Compound:inchikey:DNXHEGUUPJUMQT-CBZIJGRNSA-N': 'Conjugated estrogens', 'Compound:inchikey:VOXZDWNPVJITMN-ZBRFXRBCSA-N': 'Estradiol'}
dementia_targets = {'Disease:DOID:10652': "Alzheimer's disease", 'Disease:DOID:1307': 'Dementia', 'Disease:DOID:8725': 'Vascular dementia'}
print('=' * 80)
print('MHT → Dementia Analysis')
print('=' * 80)
print(f'\n📋 MHT Compounds ({len(mht_compounds)}):')
for cid, name in mht_compounds.items():
    if cid in G.nodes:
        deg = G.degree(cid)
        ntype = G.nodes[cid].get('node_type', '?')
        print(f'  ✓ {name:<30s}  degree={deg}  [{ntype}]')
    else:
        print(f'  ✗ {name:<30s}  NOT IN GRAPH')
print(f'\n📋 Dementia Targets ({len(dementia_targets)}):')
for did, name in dementia_targets.items():
    if did in G.nodes:
        deg = G.degree(did)
        print(f'  ✓ {name:<30s}  degree={deg}')
    else:
        print(f'  ✗ {name:<30s}  NOT IN GRAPH')
print('\n' + '=' * 80)
print('PHASE 1 — Loading dementia embedding')
print('=' * 80)
dementia_emb = pd.read_csv('${SPOKE_DATA_ROOT}/dementia_embedding.csv')
print(f'  Loaded embedding: {len(dementia_emb):,} nodes')
print(f'  Columns: {list(dementia_emb.columns)}')
emb_dict = dict(zip(dementia_emb['node'], dementia_emb['pagerank']))
print(f'\n  MHT compounds in dementia embedding:')
mht_emb_rows = []
for cid, name in mht_compounds.items():
    pr = emb_dict.get(cid, None)
    if pr is not None:
        rank = (dementia_emb['pagerank'] > pr).sum() + 1
        percentile = (1 - rank / len(dementia_emb)) * 100
        print(f'    {name:<30s}  PR={pr:.2e}  rank={rank:>6,}/{len(dementia_emb):,}  (top {100 - percentile:.2f}%)')
        mht_emb_rows.append({'compound_id': cid, 'name': name, 'pagerank': pr, 'rank': rank, 'percentile': percentile})
    else:
        print(f'    {name:<30s}  NOT IN EMBEDDING')
        mht_emb_rows.append({'compound_id': cid, 'name': name, 'pagerank': np.nan, 'rank': np.nan, 'percentile': np.nan})
df_mht_emb = pd.DataFrame(mht_emb_rows).sort_values('pagerank', ascending=False)
print('\n' + '=' * 80)
print('PHASE 2 — Gene overlap: MHT compounds ↔ Dementia targets')
print('=' * 80)

def get_gene_neighbors(G, node_id, max_hops=2):
    """BFS to collect Gene nodes within max_hops."""
    if node_id not in G.nodes:
        return set()
    genes = set()
    visited = {node_id: 0}
    queue = [(node_id, 0)]
    while queue:
        current, depth = queue.pop(0)
        if depth >= max_hops:
            continue
        for nb in G.neighbors(current):
            if nb not in visited:
                visited[nb] = depth + 1
                queue.append((nb, depth + 1))
                if nb in G.nodes and G.nodes[nb].get('node_type') == 'Gene':
                    genes.add(nb)
    return genes

def jaccard(s1, s2):
    u = s1 | s2
    return len(s1 & s2) / len(u) if u else 0.0
print(f'\n  Collecting genes for dementia targets ({MAX_HOPS} hops)...')
dementia_genes = set()
dementia_genes_per_target = {}
for did, dname in dementia_targets.items():
    genes_d = get_gene_neighbors(G, did, MAX_HOPS)
    dementia_genes_per_target[dname] = genes_d
    dementia_genes |= genes_d
    print(f'    {dname}: {len(genes_d):,} genes')
print(f'  Total dementia genes (union): {len(dementia_genes):,}')
print(f'\n  Collecting genes for MHT compounds ({MAX_HOPS} hops)...')
mht_genes = {}
mht_all_genes = set()
for cid, cname in mht_compounds.items():
    genes_c = get_gene_neighbors(G, cid, MAX_HOPS)
    mht_genes[cname] = genes_c
    mht_all_genes |= genes_c
    print(f'    {cname}: {len(genes_c):,} genes')
print(f'  Total MHT genes (union): {len(mht_all_genes):,}')
print(f'\n  Per-compound overlap with dementia genes:')
compound_overlap_rows = []
for cname, cgenes in mht_genes.items():
    overlap = cgenes & dementia_genes
    jac = jaccard(cgenes, dementia_genes)
    pct = len(overlap) / len(cgenes) * 100 if cgenes else 0
    print(f'    {cname:<30s}  {len(overlap):>4} shared / {len(cgenes):>5} total  ({pct:.1f}%)  Jaccard={jac:.4f}')
    row = {'compound': cname, 'n_genes': len(cgenes), 'overlap_all_dementia': len(overlap), 'jaccard_all_dementia': jac}
    for dname, dgenes in dementia_genes_per_target.items():
        row[f'overlap_{dname}'] = len(cgenes & dgenes)
        row[f'jaccard_{dname}'] = jaccard(cgenes, dgenes)
    compound_overlap_rows.append(row)
df_compound_overlap = pd.DataFrame(compound_overlap_rows)
mht_dementia_overlap = mht_all_genes & dementia_genes
print(f'\n  MHT (union) ∩ Dementia: {len(mht_dementia_overlap):,} genes')
print(f'  Jaccard: {jaccard(mht_all_genes, dementia_genes):.4f}')
print('\n' + '=' * 80)
print("PHASE 3 — Fisher's exact test per MHT compound")
print('=' * 80)
all_graph_genes = set((n for n in G.nodes if G.nodes[n].get('node_type') == 'Gene'))
print(f'  Background: {len(all_graph_genes):,} total genes in SPOKE\n')
fisher_rows = []
for cname, cgenes in mht_genes.items():
    a = len(cgenes & dementia_genes)
    b = len(cgenes - dementia_genes)
    c = len(dementia_genes - cgenes)
    d = len(all_graph_genes - cgenes - dementia_genes)
    table = np.array([[a, b], [c, d]])
    odds_ratio, pval = fisher_exact(table, alternative='greater')
    sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
    print(f'  {cname:<30s}  OR={odds_ratio:>6.2f}  p={pval:.2e}  {sig}')
    fisher_rows.append({'compound': cname, 'a': a, 'b': b, 'c': c, 'd': d, 'odds_ratio': odds_ratio, 'p_value': pval})
df_fisher = pd.DataFrame(fisher_rows).sort_values('odds_ratio', ascending=False)
a_u = len(mht_all_genes & dementia_genes)
b_u = len(mht_all_genes - dementia_genes)
c_u = len(dementia_genes - mht_all_genes)
d_u = len(all_graph_genes - mht_all_genes - dementia_genes)
or_u, p_u = fisher_exact(np.array([[a_u, b_u], [c_u, d_u]]), alternative='greater')
print(f"\n  MHT (union)                   OR={or_u:>6.2f}  p={p_u:.2e}  {('***' if p_u < 0.001 else '**' if p_u < 0.01 else '*' if p_u < 0.05 else 'ns')}")
print('\n' + '=' * 80)
print('PHASE 4 — Null distribution (random compounds vs dementia)')
print('=' * 80)
mht_ids = set(mht_compounds.keys())
compound_pool = [n for n in G.nodes if G.nodes[n].get('node_type') == 'Compound' and n not in mht_ids]
print(f'  Compound pool: {len(compound_pool):,} compounds (excluding MHT)')
print(f'  MHT sample size: {len(mht_compounds)}')
print(f'  Permutations: {N_PERMUTATIONS:,}')
print(f'\n  Pre-computing gene neighbors for all compounds ({MAX_HOPS} hops)...')
t0 = time.time()
compound_gene_cache = {}
for c in compound_pool:
    genes_c = set()
    if MAX_HOPS == 1:
        for nb in G.neighbors(c):
            if nb in G.nodes and G.nodes[nb].get('node_type') == 'Gene':
                genes_c.add(nb)
    else:
        visited = {c: 0}
        queue = [(c, 0)]
        while queue:
            current, depth = queue.pop(0)
            if depth >= MAX_HOPS:
                continue
            for nb in G.neighbors(current):
                if nb not in visited:
                    visited[nb] = depth + 1
                    queue.append((nb, depth + 1))
                    if nb in G.nodes and G.nodes[nb].get('node_type') == 'Gene':
                        genes_c.add(nb)
    if genes_c:
        compound_gene_cache[c] = genes_c
compound_pool_filtered = list(compound_gene_cache.keys())
t_cache = time.time() - t0
print(f'  Pre-computation done in {t_cache:.1f}s')
print(f'  Compounds with ≥1 gene: {len(compound_pool_filtered):,}')
null_overlap_counts = []
null_jaccard_vals = []
null_gene_counts = []
null_per_compound_overlaps = []
print(f'\n  Running {N_PERMUTATIONS} permutations...')
t0 = time.time()
n_mht = len(mht_compounds)
for i in range(N_PERMUTATIONS):
    sampled = random.sample(compound_pool_filtered, min(n_mht, len(compound_pool_filtered)))
    random_genes = set()
    max_single_overlap = 0
    for c in sampled:
        c_genes = compound_gene_cache[c]
        random_genes |= c_genes
        single_ovlp = len(c_genes & dementia_genes)
        max_single_overlap = max(max_single_overlap, single_ovlp)
    ovlp = len(random_genes & dementia_genes)
    jac = jaccard(random_genes, dementia_genes)
    null_overlap_counts.append(ovlp)
    null_jaccard_vals.append(jac)
    null_gene_counts.append(len(random_genes))
    null_per_compound_overlaps.append(max_single_overlap)
    if (i + 1) % 200 == 0:
        elapsed = time.time() - t0
        print(f'    [{i + 1}/{N_PERMUTATIONS}] median overlap={int(np.median(null_overlap_counts))}, ({elapsed:.1f}s)')
t_perm = time.time() - t0
print(f'\n  ✅ {N_PERMUTATIONS} permutations in {t_perm:.1f}s')
observed_overlap = len(mht_dementia_overlap)
observed_jaccard = jaccard(mht_all_genes, dementia_genes)
observed_genes = len(mht_all_genes)

def enrichment_stats(observed, null_dist, label):
    null_arr = np.array(null_dist, dtype=float)
    mean_n = np.mean(null_arr)
    std_n = np.std(null_arr)
    z = (observed - mean_n) / std_n if std_n > 0 else np.inf
    p_emp = np.mean(null_arr >= observed)
    sig = '***' if p_emp < 0.001 else '**' if p_emp < 0.01 else '*' if p_emp < 0.05 else 'ns'
    obs_str = f'{observed:,}' if isinstance(observed, int) else f'{observed:.4f}'
    print(f'  {label}:')
    print(f'    Observed     : {obs_str}')
    print(f'    Null mean±sd : {mean_n:.1f} ± {std_n:.1f}')
    print(f'    Z-score      : {z:.2f}')
    print(f'    Empirical p  : {p_emp:.4f}  {sig}')
    print()
    return {'label': label, 'observed': observed, 'null_mean': mean_n, 'null_std': std_n, 'z_score': z, 'p_empirical': p_emp}
print('\n📊 Null Distribution Results:')
print('-' * 60)
null_results = []
null_results.append(enrichment_stats(observed_overlap, null_overlap_counts, 'MHT ∩ Dementia (gene count)'))
null_results.append(enrichment_stats(observed_jaccard, null_jaccard_vals, 'Jaccard MHT–Dementia'))
null_results.append(enrichment_stats(observed_genes, null_gene_counts, 'Gene set size (sanity check)'))
print('\n' + '=' * 80)
print('PHASE 5 — Per-compound null distributions')
print('=' * 80)
per_compound_null = {}
print('  Building per-compound null distributions...')
single_overlaps_all = []
for c in compound_pool_filtered:
    ovlp = len(compound_gene_cache[c] & dementia_genes)
    single_overlaps_all.append(ovlp)
single_overlaps_arr = np.array(single_overlaps_all)
for cname, cgenes in mht_genes.items():
    obs_ovlp = len(cgenes & dementia_genes)
    p_emp = np.mean(single_overlaps_arr >= obs_ovlp)
    z = (obs_ovlp - np.mean(single_overlaps_arr)) / np.std(single_overlaps_arr) if np.std(single_overlaps_arr) > 0 else 0
    sig = '***' if p_emp < 0.001 else '**' if p_emp < 0.01 else '*' if p_emp < 0.05 else 'ns'
    print(f'  {cname:<30s}  overlap={obs_ovlp:>4}  null_mean={np.mean(single_overlaps_arr):.1f}  z={z:.2f}  p={p_emp:.4f}  {sig}')
    per_compound_null[cname] = {'observed': obs_ovlp, 'z': z, 'p': p_emp, 'null_mean': np.mean(single_overlaps_arr), 'null_std': np.std(single_overlaps_arr)}
print('\n  Embedding rank comparison:')
compound_prs = []
for c in compound_pool_filtered:
    pr = emb_dict.get(c, None)
    if pr is not None:
        compound_prs.append(pr)
compound_prs_arr = np.array(compound_prs)
for cid, cname in mht_compounds.items():
    pr = emb_dict.get(cid, None)
    if pr is not None:
        p_emp_pr = np.mean(compound_prs_arr >= pr)
        z_pr = (pr - np.mean(compound_prs_arr)) / np.std(compound_prs_arr) if np.std(compound_prs_arr) > 0 else 0
        sig = '***' if p_emp_pr < 0.001 else '**' if p_emp_pr < 0.01 else '*' if p_emp_pr < 0.05 else 'ns'
        print(f'  {cname:<30s}  PR={pr:.2e}  z={z_pr:.2f}  p={p_emp_pr:.4f}  {sig}')
print('\n' + '=' * 80)
print('PHASE 6 — 🔒 Saving results')
print('=' * 80)
suffix = f'_{MAX_HOPS}hops'
df_mht_emb.to_csv(f'{SAVE_PATH}mht_dementia_embedding_ranks{suffix}.csv', index=False)
print(f'  💾 mht_dementia_embedding_ranks{suffix}.csv')
df_compound_overlap.to_csv(f'{SAVE_PATH}mht_dementia_compound_overlap{suffix}.csv', index=False)
print(f'  💾 mht_dementia_compound_overlap{suffix}.csv')
df_fisher.to_csv(f'{SAVE_PATH}mht_dementia_fisher{suffix}.csv', index=False)
print(f'  💾 mht_dementia_fisher{suffix}.csv')
df_null = pd.DataFrame(null_results)
df_null.to_csv(f'{SAVE_PATH}mht_dementia_null_results{suffix}.csv', index=False)
print(f'  💾 mht_dementia_null_results{suffix}.csv')
df_per_comp_null = pd.DataFrame([{'compound': k, **v} for k, v in per_compound_null.items()])
df_per_comp_null.to_csv(f'{SAVE_PATH}mht_dementia_per_compound_null{suffix}.csv', index=False)
print(f'  💾 mht_dementia_per_compound_null{suffix}.csv')
mht_dem_genes_rows = []
for gid in sorted(mht_dementia_overlap):
    gname = G.nodes[gid].get('name', gid) if gid in G.nodes else gid
    sources = []
    for cname, cgenes in mht_genes.items():
        if gid in cgenes:
            sources.append(cname)
    mht_dem_genes_rows.append({'gene_id': gid, 'gene_name': gname, 'mht_sources': '; '.join(sources), 'n_mht_sources': len(sources)})
df_shared_genes = pd.DataFrame(mht_dem_genes_rows)
df_shared_genes.to_csv(f'{SAVE_PATH}mht_dementia_shared_genes{suffix}.csv', index=False)
print(f'  💾 mht_dementia_shared_genes{suffix}.csv ({len(df_shared_genes)} genes)')
print(f'\n✅ All results saved to {SAVE_PATH}')
print('\n' + '=' * 80)
print('PHASE 7 — Generating figures')
print('=' * 80)
fig1, (ax1a, ax1b) = plt.subplots(1, 2, figsize=(16, 5))
heat_cols = [c for c in df_compound_overlap.columns if c.startswith('overlap_')]
heat_data = df_compound_overlap.set_index('compound')[heat_cols].copy()
heat_data.columns = [c.replace('overlap_', '') for c in heat_cols]
sns.heatmap(heat_data, annot=True, fmt='d', cmap='OrRd', ax=ax1a, linewidths=0.5, linecolor='white', cbar_kws={'label': 'Shared genes'})
ax1a.set_title('A. Gene Overlap with Dementia Targets', fontsize=13, fontweight='bold')
ax1a.set_ylabel('')
jac_cols = [c for c in df_compound_overlap.columns if c.startswith('jaccard_')]
jac_data = df_compound_overlap.set_index('compound')[jac_cols].copy()
jac_data.columns = [c.replace('jaccard_', '') for c in jac_cols]
sns.heatmap(jac_data, annot=True, fmt='.4f', cmap='YlGnBu', ax=ax1b, linewidths=0.5, linecolor='white', cbar_kws={'label': 'Jaccard index'})
ax1b.set_title('B. Jaccard Similarity', fontsize=13, fontweight='bold')
ax1b.set_ylabel('')
fig1.suptitle(f'MHT Compound–Dementia Gene Overlap ({MAX_HOPS} hops)', fontsize=14, fontweight='bold')
plt.tight_layout()
fig1.savefig(f'{SAVE_PATH}mht_dementia_overlap_heatmap{suffix}.png', dpi=300, bbox_inches='tight')
fig1.savefig(f'{SAVE_PATH}mht_dementia_overlap_heatmap{suffix}.pdf', dpi=300, bbox_inches='tight')
print('  ✅ Overlap heatmap saved')
plt.show()
fig2, axes2 = plt.subplots(1, 3, figsize=(18, 5))

def plot_null(ax, null_dist, observed, label, color):
    null_arr = np.array(null_dist, dtype=float)
    p_emp = np.mean(null_arr >= observed)
    z = (observed - np.mean(null_arr)) / np.std(null_arr) if np.std(null_arr) > 0 else 0
    ax.hist(null_arr, bins=50, color=COLORS['random'], edgecolor='gray', alpha=0.7, density=True, label='Random compounds')
    tail = null_arr[null_arr >= observed]
    if len(tail) > 0:
        ax.hist(tail, bins=50, color=color, alpha=0.3, density=True)
    obs_label = f'{observed:,}' if isinstance(observed, int) else f'{observed:.4f}'
    ax.axvline(observed, color=color, linewidth=2.5, linestyle='--', label=f'MHT = {obs_label}')
    sig = '***' if p_emp < 0.001 else '**' if p_emp < 0.01 else '*' if p_emp < 0.05 else 'ns'
    ax.text(0.97, 0.95, f'p = {p_emp:.4f} {sig}\nz = {z:.1f}', transform=ax.transAxes, ha='right', va='top', fontsize=10, bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    ax.set_xlabel(label, fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.legend(fontsize=9, loc='upper left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
plot_null(axes2[0], null_overlap_counts, observed_overlap, 'Shared genes with Dementia', COLORS['mht'])
axes2[0].set_title('A. MHT ∩ Dementia (gene count)', fontsize=13, fontweight='bold')
plot_null(axes2[1], null_jaccard_vals, observed_jaccard, 'Jaccard index', COLORS['mht'])
axes2[1].set_title('B. Jaccard MHT–Dementia', fontsize=13, fontweight='bold')
plot_null(axes2[2], null_gene_counts, observed_genes, 'Gene set size', COLORS['mht'])
axes2[2].set_title('C. Gene set size (sanity)', fontsize=13, fontweight='bold')
fig2.suptitle(f'Null Distribution: MHT vs {N_PERMUTATIONS:,} Random Compound Sets\n({MAX_HOPS} hops, {len(compound_pool_filtered):,} compounds in pool)', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
fig2.savefig(f'{SAVE_PATH}mht_dementia_null_distribution{suffix}.png', dpi=300, bbox_inches='tight')
fig2.savefig(f'{SAVE_PATH}mht_dementia_null_distribution{suffix}.pdf', dpi=300, bbox_inches='tight')
print('  ✅ Null distribution figure saved')
plt.show()
fig3, ax3 = plt.subplots(figsize=(12, 6))
compound_names = list(mht_genes.keys())
compound_overlaps = [len(mht_genes[c] & dementia_genes) for c in compound_names]
compound_colors = [COLORS.get(c.split()[0].lower(), COLORS['mht']) for c in compound_names]
color_map = {'Conjugated estrogens': COLORS['conj_estrogens'], 'Estradiol': COLORS['estradiol'], 'Progesterone': COLORS['progesterone'], 'Medroxyprogesterone acetate': COLORS['medroxy'], 'Norethindrone acetate': COLORS['norethindrone']}
bar_colors = [color_map.get(c, COLORS['mht']) for c in compound_names]
bars = ax3.bar(range(len(compound_names)), compound_overlaps, color=bar_colors, edgecolor='white', linewidth=2, alpha=0.85)
null_mean = np.mean(single_overlaps_all)
null_p95 = np.percentile(single_overlaps_all, 95)
null_p99 = np.percentile(single_overlaps_all, 99)
ax3.axhline(null_mean, color='gray', linestyle='-', linewidth=1.5, label=f'Random mean = {null_mean:.0f}')
ax3.axhline(null_p95, color='orange', linestyle='--', linewidth=1, label=f'95th percentile = {null_p95:.0f}')
ax3.axhline(null_p99, color='red', linestyle='--', linewidth=1, label=f'99th percentile = {null_p99:.0f}')
for bar, ovlp, cname in zip(bars, compound_overlaps, compound_names):
    h = bar.get_height()
    p = per_compound_null[cname]['p']
    sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    ax3.text(bar.get_x() + bar.get_width() / 2.0, h + max(compound_overlaps) * 0.02, f'{ovlp}\n{sig}', ha='center', va='bottom', fontweight='bold', fontsize=10)
ax3.set_xticks(range(len(compound_names)))
ax3.set_xticklabels(compound_names, rotation=25, ha='right', fontsize=10)
ax3.set_ylabel('Genes shared with Dementia', fontsize=12, fontweight='bold')
ax3.set_title('Per-MHT Compound: Dementia Gene Overlap vs Random Compound Null', fontsize=14, fontweight='bold')
ax3.legend(loc='upper right', fontsize=9)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.set_ylim(0, max(max(compound_overlaps), null_p99) * 1.3)
plt.tight_layout()
fig3.savefig(f'{SAVE_PATH}mht_dementia_per_compound{suffix}.png', dpi=300, bbox_inches='tight')
fig3.savefig(f'{SAVE_PATH}mht_dementia_per_compound{suffix}.pdf', dpi=300, bbox_inches='tight')
print('  ✅ Per-compound bar chart saved')
plt.show()
fig4, ax4 = plt.subplots(figsize=(12, 6))
df_fisher_sorted = df_fisher.sort_values('odds_ratio', ascending=False)
bar_colors_f = [color_map.get(c, COLORS['mht']) for c in df_fisher_sorted['compound']]
bars_f = ax4.bar(range(len(df_fisher_sorted)), df_fisher_sorted['odds_ratio'], color=bar_colors_f, edgecolor='white', linewidth=2, alpha=0.85)
ax4.axhline(1.0, color='gray', linestyle='--', linewidth=1.5, label='OR = 1 (no enrichment)')
for bar, (_, row) in zip(bars_f, df_fisher_sorted.iterrows()):
    h = bar.get_height()
    sig = '***' if row['p_value'] < 0.001 else '**' if row['p_value'] < 0.01 else '*' if row['p_value'] < 0.05 else 'ns'
    ax4.text(bar.get_x() + bar.get_width() / 2.0, h + max(df_fisher_sorted['odds_ratio']) * 0.02, f"OR={row['odds_ratio']:.1f}\np={row['p_value']:.1e}\n{sig}", ha='center', va='bottom', fontsize=9, fontweight='bold')
ax4.set_xticks(range(len(df_fisher_sorted)))
ax4.set_xticklabels(df_fisher_sorted['compound'], rotation=25, ha='right', fontsize=10)
ax4.set_ylabel('Odds Ratio', fontsize=12, fontweight='bold')
ax4.set_title("Fisher's Exact Test: MHT Compound Gene Enrichment in Dementia", fontsize=14, fontweight='bold')
ax4.legend(fontsize=10)
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
plt.tight_layout()
fig4.savefig(f'{SAVE_PATH}mht_dementia_fisher{suffix}.png', dpi=300, bbox_inches='tight')
fig4.savefig(f'{SAVE_PATH}mht_dementia_fisher{suffix}.pdf', dpi=300, bbox_inches='tight')
print("  ✅ Fisher's exact test figure saved")
plt.show()
fig5, ax5 = plt.subplots(figsize=(10, 6))
ax5.hist(compound_prs_arr, bins=100, color=COLORS['random'], edgecolor='gray', alpha=0.6, density=True, label=f'Random compounds (n={len(compound_prs_arr):,})')
for cid, cname in mht_compounds.items():
    pr = emb_dict.get(cid, None)
    if pr is not None:
        ax5.axvline(pr, color=color_map.get(cname, COLORS['mht']), linewidth=2, linestyle='--', label=f'{cname} ({pr:.2e})')
ax5.set_xlabel('Dementia PageRank Score', fontsize=12, fontweight='bold')
ax5.set_ylabel('Density', fontsize=12)
ax5.set_title('MHT Compound Proximity to Dementia (PPR Embedding)', fontsize=14, fontweight='bold')
ax5.legend(fontsize=8, loc='upper right')
ax5.spines['top'].set_visible(False)
ax5.spines['right'].set_visible(False)
plt.tight_layout()
fig5.savefig(f'{SAVE_PATH}mht_dementia_embedding_ranks{suffix}.png', dpi=300, bbox_inches='tight')
fig5.savefig(f'{SAVE_PATH}mht_dementia_embedding_ranks{suffix}.pdf', dpi=300, bbox_inches='tight')
print('  ✅ Embedding rank figure saved')
plt.show()
from matplotlib.gridspec import GridSpec
fig6 = plt.figure(figsize=(20, 12))
gs6 = GridSpec(2, 3, figure=fig6, hspace=0.35, wspace=0.35)
ax6a = fig6.add_subplot(gs6[0, 0])
bars_a = ax6a.barh(range(len(compound_names)), compound_overlaps, color=[color_map.get(c, COLORS['mht']) for c in compound_names], edgecolor='white', alpha=0.85)
ax6a.axvline(null_mean, color='gray', linestyle='--', linewidth=1, label='Null mean')
ax6a.axvline(null_p95, color='red', linestyle='--', linewidth=1, label='95th %ile')
ax6a.set_yticks(range(len(compound_names)))
ax6a.set_yticklabels(compound_names, fontsize=9)
ax6a.set_xlabel('Shared genes with Dementia', fontsize=10)
ax6a.set_title('A. Per-compound Overlap', fontsize=12, fontweight='bold')
ax6a.legend(fontsize=8)
ax6a.spines['top'].set_visible(False)
ax6a.spines['right'].set_visible(False)
ax6b = fig6.add_subplot(gs6[0, 1])
or_vals = [df_fisher[df_fisher['compound'] == c]['odds_ratio'].values[0] if c in df_fisher['compound'].values else 0 for c in compound_names]
bars_b = ax6b.barh(range(len(compound_names)), or_vals, color=[color_map.get(c, COLORS['mht']) for c in compound_names], edgecolor='white', alpha=0.85)
ax6b.axvline(1.0, color='gray', linestyle='--', linewidth=1.5)
ax6b.set_yticks(range(len(compound_names)))
ax6b.set_yticklabels(compound_names, fontsize=9)
ax6b.set_xlabel('Odds Ratio', fontsize=10)
ax6b.set_title("B. Fisher's Exact Test (OR)", fontsize=12, fontweight='bold')
ax6b.spines['top'].set_visible(False)
ax6b.spines['right'].set_visible(False)
ax6c = fig6.add_subplot(gs6[0, 2])
pr_vals = [emb_dict.get(cid, 0) for cid in mht_compounds.keys()]
pr_names = list(mht_compounds.values())
bars_c = ax6c.barh(range(len(pr_names)), pr_vals, color=[color_map.get(c, COLORS['mht']) for c in pr_names], edgecolor='white', alpha=0.85)
ax6c.axvline(np.median(compound_prs_arr), color='gray', linestyle='--', linewidth=1, label='Median compound')
ax6c.set_yticks(range(len(pr_names)))
ax6c.set_yticklabels(pr_names, fontsize=9)
ax6c.set_xlabel('Dementia PageRank', fontsize=10)
ax6c.set_title('C. Proximity to Dementia (PPR)', fontsize=12, fontweight='bold')
ax6c.legend(fontsize=8)
ax6c.spines['top'].set_visible(False)
ax6c.spines['right'].set_visible(False)
ax6d = fig6.add_subplot(gs6[1, 0])
plot_null(ax6d, null_overlap_counts, observed_overlap, 'Shared genes', COLORS['mht'])
ax6d.set_title('D. MHT vs Random (overlap)', fontsize=12, fontweight='bold')
ax6e = fig6.add_subplot(gs6[1, 1])
plot_null(ax6e, null_jaccard_vals, observed_jaccard, 'Jaccard', COLORS['mht'])
ax6e.set_title('E. MHT vs Random (Jaccard)', fontsize=12, fontweight='bold')
ax6f = fig6.add_subplot(gs6[1, 2])
plot_null(ax6f, null_gene_counts, observed_genes, 'Gene set size', COLORS['mht'])
ax6f.set_title('F. Gene set size (sanity)', fontsize=12, fontweight='bold')
fig6.suptitle('MHT → Dementia: Comprehensive Analysis Summary', fontsize=16, fontweight='bold', y=1.01)
plt.tight_layout()
fig6.savefig(f'{SAVE_PATH}mht_dementia_summary{suffix}.png', dpi=300, bbox_inches='tight')
fig6.savefig(f'{SAVE_PATH}mht_dementia_summary{suffix}.pdf', dpi=300, bbox_inches='tight')
print('  ✅ Summary figure saved')
plt.show()
if USE_DWPC:
    print('\n' + '=' * 80)
    print('PHASE 8 — DWPC analysis (MHT compounds → Dementia)')
    print('=' * 80)
    print('  This may take a while...\n')
    try:
        dwpc_results = []
        for cid, cname in mht_compounds.items():
            if cid not in G.nodes:
                print(f'  ⚠️  {cname} not in graph, skipping')
                continue
            source_type = G.nodes[cid].get('node_type')
            for did, dname in dementia_targets.items():
                if did not in G.nodes:
                    continue
                target_type = G.nodes[did].get('node_type')
                try:
                    metapaths = calculator.find_metapaths(source_type, target_type, max_length=DWPC_MAX_LENGTH)
                    for metapath in metapaths:
                        dwpc = calculator.compute_dwpc(cid, did, metapath)
                        path_count = calculator.compute_simple_path_count(cid, did, metapath)
                        dwpc_results.append({'source_id': cid, 'source_name': cname, 'target_id': did, 'target_name': dname, 'metapath': ' → '.join(metapath), 'length': len(metapath), 'dwpc': dwpc, 'path_count': path_count})
                except Exception as e:
                    print(f'  ❌ Error {cname}→{dname}: {e}')
                    continue
            print(f"  ✓ {cname}: {sum((1 for r in dwpc_results if r['source_name'] == cname))} metapaths")
        df_dwpc = pd.DataFrame(dwpc_results)
        if len(df_dwpc) > 0:
            df_dwpc = df_dwpc.sort_values('dwpc', ascending=False).reset_index(drop=True)
            df_dwpc.to_csv(f'{SAVE_PATH}mht_dementia_dwpc{suffix}.csv', index=False)
            print(f'\n  💾 mht_dementia_dwpc{suffix}.csv ({len(df_dwpc)} rows)')
            print(f'\n  Top 20 DWPC results:')
            for _, row in df_dwpc.head(20).iterrows():
                print(f"    {row['source_name']:<25s} → {row['target_name']:<20s}  DWPC={row['dwpc']:.4f}  paths={row['path_count']}  ({row['metapath']})")
            best_dwpc = df_dwpc.groupby(['source_name', 'target_name'])['dwpc'].max().reset_index()
            dwpc_pivot = best_dwpc.pivot(index='source_name', columns='target_name', values='dwpc')
            fig_dwpc, ax_dwpc = plt.subplots(figsize=(10, 5))
            sns.heatmap(dwpc_pivot, annot=True, fmt='.4f', cmap='Purples', ax=ax_dwpc, linewidths=0.5, linecolor='white', cbar_kws={'label': 'Max DWPC score'})
            ax_dwpc.set_title('MHT → Dementia: Best DWPC Score per Compound–Target Pair', fontsize=13, fontweight='bold')
            ax_dwpc.set_ylabel('')
            ax_dwpc.set_xlabel('')
            plt.tight_layout()
            fig_dwpc.savefig(f'{SAVE_PATH}mht_dementia_dwpc_heatmap{suffix}.png', dpi=300, bbox_inches='tight')
            fig_dwpc.savefig(f'{SAVE_PATH}mht_dementia_dwpc_heatmap{suffix}.pdf', dpi=300, bbox_inches='tight')
            print('  ✅ DWPC heatmap saved')
            plt.show()
        else:
            print('  ⚠️  No DWPC results found')
    except NameError:
        print('  ⚠️  NetworkDWPC calculator not defined. Skipping DWPC analysis.')
        print('     Create it with: calculator = NetworkDWPC(G)')
    except Exception as e:
        print(f'  ❌ DWPC failed: {e}')
        print('     Core results are already saved.')
print('\n' + '=' * 80)
print('🏁 MHT → DEMENTIA ANALYSIS COMPLETE')
print('=' * 80)
print(f"\nConfiguration:\n  • {MAX_HOPS} hops, {N_PERMUTATIONS} permutations\n  • {len(mht_compounds)} MHT compounds, {len(dementia_targets)} dementia targets\n  • {len(compound_pool_filtered):,} random compounds in pool\n\nResults:\n  • MHT union genes: {len(mht_all_genes):,}\n  • Dementia genes: {len(dementia_genes):,}\n  • MHT ∩ Dementia: {len(mht_dementia_overlap):,}\n\nFisher's (MHT union → Dementia): OR={or_u:.2f}, p={p_u:.2e}\n\nPer-compound highlights:")
for cname in compound_names:
    pcn = per_compound_null[cname]
    fr = df_fisher[df_fisher['compound'] == cname].iloc[0]
    ovlp = len(mht_genes[cname] & dementia_genes)
    print(f"  • {cname:<30s}  genes={len(mht_genes[cname]):>5}  ∩dementia={ovlp:>4}  OR={fr['odds_ratio']:.1f}  p_fisher={fr['p_value']:.1e}  p_null={pcn['p']:.4f}")
print(f'\n  All files saved to: {SAVE_PATH}')
