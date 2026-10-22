# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""
Recompute gene_summary for PM (DOID:10787) → Dementia (DOID:1307) with MAX_PATHS=50
restricted to the first TOP_METAPATHS_LIMIT metapaths that have ≥1 path in the real analysis.
Saves pm_dementia_dwpc_recurring_genes_v2.csv.
"""
import os
import os
from collections import Counter, defaultdict
from typing import List
import numpy as np
import pandas as pd
try:
    import networkx as nx
except Exception:
    nx = None
SOURCE_NODE = 'Disease:DOID:10787'
TARGET_NODE = 'Disease:DOID:1307'
MAX_PATHS = int(os.environ.get('MAX_PATHS', '50'))
TOP_METAPATHS_LIMIT = int(os.environ.get('TOP_METAPATHS_LIMIT', '11'))
BASE = '${SPOKE_DATA_ROOT}'
GFILE = os.environ.get('GFILE', os.path.join(BASE, 'spoke_graph.gpickle'))
METAPATHS_CSV = os.path.join(BASE, 'df_premature_menopause_to_dementia.csv')
OUTCSV = os.path.join(BASE, 'pm_dementia_dwpc_recurring_genes_v2.csv')

def extract_intermediates_from_paths(G, source_node, target_node, metapath: List[str], max_paths: int=50):
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
            for _, edata in G.get_edge_data(node, nbr).items():
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
                inter_nodes.append(n)
    return (inter_nodes, len(results))

def main():
    if nx is None:
        raise SystemExit('networkx not available on PATH')
    if not os.path.exists(GFILE):
        raise SystemExit(f'Graph not found. Please save G to {GFILE} from the notebook.')
    if not os.path.exists(METAPATHS_CSV):
        raise SystemExit(f'Missing metapaths CSV: {METAPATHS_CSV}')
    try:
        G = nx.read_gpickle(GFILE) if hasattr(nx, 'read_gpickle') else None
    except Exception:
        G = None
    if G is None:
        import pickle
        with open(GFILE, 'rb') as fh:
            G = pickle.load(fh)
    df_mp = pd.read_csv(METAPATHS_CSV)
    df_mp = df_mp[df_mp.get('path_count', 0) > 0].sort_values('dwpc', ascending=False).head(TOP_METAPATHS_LIMIT)
    gene_occ = Counter()
    gene_metapaths = defaultdict(set)
    for _, row in df_mp.iterrows():
        schema = str(row['metapath']).split(' → ')
        inters, n_paths = extract_intermediates_from_paths(G, SOURCE_NODE, TARGET_NODE, schema, max_paths=MAX_PATHS)
        for g in inters:
            gene_occ[g] += 1
            gene_metapaths[g].add(row['metapath'])
    if not gene_occ:
        raise SystemExit('No genes recovered with current settings (check metapaths or MAX_PATHS)')
    out = pd.DataFrame({'gene_name': list(gene_occ.keys()), 'total_occurrences': [gene_occ[g] for g in gene_occ.keys()], 'n_metapaths': [len(gene_metapaths[g]) for g in gene_occ.keys()]}).sort_values(['total_occurrences', 'n_metapaths', 'gene_name'], ascending=[False, False, True])
    out.to_csv(OUTCSV, index=False)
    print(f'Saved {OUTCSV} (n={len(out)})')
if __name__ == '__main__':
    main()
