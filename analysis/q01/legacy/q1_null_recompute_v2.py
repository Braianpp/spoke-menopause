# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""
Recompute concept-swap null for Premature Menopause (DOID:10787) → Dementia (DOID:1307)
with symmetric path cap and restricted to metapaths that have ≥1 path in the real analysis.

Assumptions / prerequisites (same as notebook):
- G: networkx MultiGraph (undirected accepted) with edge attribute 'predicate'
- gene_summary: DataFrame with columns ['gene_name','total_occurrences','n_metapaths']
- df_premature_menopause_to_dementia (or similar): has columns ['metapath','path_count','dwpc']
- If extract_intermediates_from_paths is not defined in the session, a fallback DFS is provided here.

Usage options:
- Import and run inside the dwpc.ipynb session where G, gene_summary, and df_* already exist:


- Or run as a script if you also provide a pickled G and CSVs for df and genes.
"""
import os
import os
import math
import random
from collections import Counter
from typing import List, Tuple
import numpy as np
import pandas as pd
try:
    import networkx as nx
except Exception:
    nx = None
SOURCE_NODE = 'Disease:DOID:10787'
TARGET_NODE = 'Disease:DOID:1307'
N_RANDOM_PAIRS = int(os.environ.get('N_RANDOM_PAIRS', '1500'))
MAX_PATHS = int(os.environ.get('MAX_PATHS', '50'))
TOP_METAPATHS_LIMIT = int(os.environ.get('TOP_METAPATHS_LIMIT', '11'))
ALPHA = float(os.environ.get('ALPHA', '0.05'))
SEED = int(os.environ.get('SEED', '42'))
OUTDIR = os.environ.get('OUTDIR', '${SPOKE_DATA_ROOT}/biorouter/figs')
TABLE_DIR = os.environ.get('TABLE_DIR', '${SPOKE_DATA_ROOT}')
NULL_TABLE = os.path.join(TABLE_DIR, 'pm_dementia_null_gene_pvalues_v2.csv')
YIELD_TABLE = os.path.join(TABLE_DIR, 'null_distribution_results_v2.csv')
random.seed(SEED)
np.random.seed(SEED)

def extract_intermediates_from_paths(G, source_node, target_node, metapath: List[str], max_paths: int=50):
    """DFS walk constrained by a sequence of edge predicates (undirected ok).
    Returns a list of dicts: [{'node_type': 'Gene', 'node_name': ...}, ...]
    And the number of paths explored (<= max_paths).
    """
    results = []
    path = [source_node]

    def step(node, i):
        if i == len(metapath):
            if node == target_node:
                results.append(path.copy())
            return
        if len(results) >= max_paths:
            return
        required_pred = metapath[i]
        for nbr in G.neighbors(node):
            for key, edata in G.get_edge_data(node, nbr).items():
                if edata.get('predicate') == required_pred:
                    path.append(nbr)
                    step(nbr, i + 1)
                    path.pop()
                    if len(results) >= max_paths:
                        return
    step(source_node, 0)
    inter_nodes = []
    for p in results:
        for n in p[1:-1]:
            nd = G.nodes[n]
            if nd.get('node_type') == 'Gene':
                inter_nodes.append({'node_type': 'Gene', 'node_name': n})
    return (inter_nodes, len(results))

def get_node_predicates(G, node) -> set:
    preds = set()
    for nbr in G.neighbors(node):
        for _, edata in G.get_edge_data(node, nbr).items():
            pred = edata.get('predicate')
            if pred:
                preds.add(pred)
    return preds

def run_null(G=None, gene_summary: pd.DataFrame=None, df_metapaths: pd.DataFrame=None):
    if G is None or gene_summary is None or df_metapaths is None:
        raise RuntimeError('G, gene_summary, and df_metapaths must be provided in the session or passed explicitly.')
    df_mp = df_metapaths.copy().loc[df_metapaths['path_count'] > 0].sort_values('dwpc', ascending=False).head(TOP_METAPATHS_LIMIT)
    metapath_schemas = [mp.split(' → ') for mp in df_mp['metapath'].tolist()]
    required_source_preds = {schema[0] for schema in metapath_schemas}
    required_target_preds = {schema[-1] for schema in metapath_schemas}
    disease_nodes = [n for n, d in G.nodes(data=True) if str(d.get('node_type', '')).lower() == 'disease']
    pool = [n for n in disease_nodes if n not in (SOURCE_NODE, TARGET_NODE)]
    pred_map = {n: get_node_predicates(G, n) for n in pool}
    source_candidates = [n for n in pool if pred_map[n] & required_source_preds]
    target_candidates = [n for n in pool if pred_map[n] & required_target_preds]
    all_pairs = [(s, t) for s in source_candidates for t in target_candidates if s != t]
    random.shuffle(all_pairs)
    random_pairs = all_pairs[:min(N_RANDOM_PAIRS, len(all_pairs))]
    if len(random_pairs) < N_RANDOM_PAIRS:
        while len(random_pairs) < N_RANDOM_PAIRS:
            s = random.choice(source_candidates)
            t = random.choice(target_candidates)
            if s == t:
                continue
            random_pairs.append((s, t))
    background_gene_counts: List[Counter] = []
    for rs, rt in random_pairs:
        ctr = Counter()
        for schema in metapath_schemas:
            try:
                inter_list, n_paths = extract_intermediates_from_paths(G, rs, rt, schema, max_paths=MAX_PATHS)
            except Exception:
                continue
            for item in inter_list:
                if item.get('node_type') == 'Gene':
                    ctr[item['node_name']] += 1
        background_gene_counts.append(ctr)
    real_gene_counts = dict(zip(gene_summary['gene_name'], gene_summary['total_occurrences']))
    null_results = []
    for gene_name, real_count in real_gene_counts.items():
        bg_counts = [bg.get(gene_name, 0) for bg in background_gene_counts]
        bg_mean = float(np.mean(bg_counts))
        bg_std = float(np.std(bg_counts))
        bg_max = int(np.max(bg_counts)) if len(bg_counts) else 0
        n_geq = sum((1 for c in bg_counts if c >= real_count))
        empirical_pval = (n_geq + 1) / (len(background_gene_counts) + 1)
        z_score = (real_count - bg_mean) / bg_std if bg_std > 0 else 0.0 if real_count == bg_mean else float('inf')
        n_pairs_present = sum((1 for c in bg_counts if c > 0))
        bg_prevalence = n_pairs_present / max(1, len(background_gene_counts))
        n_metapaths = int(gene_summary.loc[gene_summary['gene_name'] == gene_name, 'n_metapaths'].values[0])
        null_results.append({'gene_name': gene_name, 'real_count': int(real_count), 'bg_mean': round(bg_mean, 2), 'bg_std': round(bg_std, 2), 'bg_max': bg_max, 'bg_prevalence': round(bg_prevalence, 3), 'z_score': round(z_score, 2) if np.isfinite(z_score) else float('inf'), 'empirical_pval': round(empirical_pval, 6), 'n_metapaths': n_metapaths})
    null_df = pd.DataFrame(null_results).sort_values('empirical_pval').reset_index(drop=True)
    m = len(null_df)
    null_df['bh_rank'] = np.arange(1, m + 1)
    null_df['bh_threshold'] = null_df['bh_rank'] / m * ALPHA
    null_df['significant_bh'] = null_df['empirical_pval'] <= null_df['bh_threshold']
    os.makedirs(TABLE_DIR, exist_ok=True)
    null_df.to_csv(NULL_TABLE, index=False)
    bg_genes_per_pair = [len(bg) for bg in background_gene_counts]
    pd.DataFrame({'bg_unique_genes': bg_genes_per_pair}).to_csv(YIELD_TABLE, index=False)
    print(f'Saved null gene table: {NULL_TABLE}')
    print(f'Saved background yield table: {YIELD_TABLE}')
    return null_df
if __name__ == '__main__':
    print(f'This script is intended to be run within the notebook session where G, gene_summary, and df_metapaths are defined.\nImport and call run_null(G, gene_summary, df_metapaths).\nCurrent settings → N_RANDOM_PAIRS={N_RANDOM_PAIRS}, MAX_PATHS={MAX_PATHS}, TOP_METAPATHS_LIMIT={TOP_METAPATHS_LIMIT}, ALPHA={ALPHA}')
