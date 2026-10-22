# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
RANDOM_SEED = 42
N_PERM = 1000
N_BINS = 20
rng = np.random.default_rng(RANDOM_SEED)
all_graph_genes = sorted([n for n in G.nodes if G.nodes[n].get('node_type') == 'Gene'])
connected_gene_universe = sorted(set(connected_genes_dict['APO'].keys()) | set(connected_genes_dict['Menopause'].keys()) | set(connected_genes_dict['Cognitive'].keys()))
USE_CONNECTED_UNIVERSE = True
gene_universe = connected_gene_universe if USE_CONNECTED_UNIVERSE else all_graph_genes
print(f"Using background universe: {('connected 2-hop gene universe' if USE_CONNECTED_UNIVERSE else 'all SPOKE genes')}")
print(f'Background genes: {len(gene_universe):,}')
degree_df = pd.DataFrame({'gene_id': gene_universe, 'degree': [G.degree(g) for g in gene_universe]})
degree_df['log_degree'] = np.log1p(degree_df['degree'])
degree_df['degree_bin'] = pd.qcut(degree_df['log_degree'], q=N_BINS, labels=False, duplicates='drop')
degree_bin_lookup = dict(zip(degree_df['gene_id'], degree_df['degree_bin']))
bin_to_genes = {b: set(sub['gene_id']) for b, sub in degree_df.groupby('degree_bin')}
print('\nDegree bins:')
print(degree_df['degree_bin'].value_counts().sort_index().to_string())

def degree_matched_sample(observed_genes, exclude=None):
    """
    Sample a random gene set with approximately the same degree-bin distribution
    as observed_genes.
    """
    if exclude is None:
        exclude = set()
    else:
        exclude = set(exclude)
    observed_genes = [g for g in observed_genes if g in degree_bin_lookup]
    obs_bins = pd.Series([degree_bin_lookup[g] for g in observed_genes]).value_counts()
    sampled = []
    for b, n_needed in obs_bins.items():
        candidates = list(bin_to_genes[b] - exclude)
        replace = len(candidates) < n_needed
        if len(candidates) == 0:
            continue
        sampled.extend(rng.choice(candidates, size=n_needed, replace=replace).tolist())
    return set(sampled)

def overlap_metrics(a, m, c):
    """
    Compute the same overlap statistics used in the main analysis.
    """
    return {'apo_meno': len(a & m), 'apo_cog': len(a & c), 'meno_cog': len(m & c), 'all_three': len(a & m & c), 'jaccard_apo_meno': jaccard(a, m), 'jaccard_apo_cog': jaccard(a, c), 'jaccard_meno_cog': jaccard(m, c)}
observed_metrics = overlap_metrics(apo_genes, meno_genes, cog_genes)
print('\nObserved overlap metrics:')
for k, v in observed_metrics.items():
    if 'jaccard' in k:
        print(f'  {k}: {v:.4f}')
    else:
        print(f'  {k}: {v:,}')
null_rows = []
for i in range(N_PERM):
    rand_apo = degree_matched_sample(apo_genes)
    rand_meno = degree_matched_sample(meno_genes)
    rand_cog = degree_matched_sample(cog_genes)
    row = overlap_metrics(rand_apo, rand_meno, rand_cog)
    row['perm'] = i + 1
    null_rows.append(row)
df_null = pd.DataFrame(null_rows)
summary_rows = []
for metric, obs in observed_metrics.items():
    null_vals = df_null[metric].values
    emp_p = (np.sum(null_vals >= obs) + 1) / (len(null_vals) + 1)
    null_mean = np.mean(null_vals)
    null_sd = np.std(null_vals, ddof=1)
    z = (obs - null_mean) / null_sd if null_sd > 0 else np.nan
    summary_rows.append({'metric': metric, 'observed': obs, 'null_mean': null_mean, 'null_sd': null_sd, 'z_score': z, 'empirical_p': emp_p})
df_null_summary = pd.DataFrame(summary_rows)
print('\nDegree-matched permutation null summary:')
display(df_null_summary)
df_null.to_csv(SAVE_PATH / f'degree_matched_null_permutations{suffix}.csv', index=False)
df_null_summary.to_csv(SAVE_PATH / f'degree_matched_null_summary{suffix}.csv', index=False)
print(f'\nSaved:')
print(f'  degree_matched_null_permutations{suffix}.csv')
print(f'  degree_matched_null_summary{suffix}.csv')
