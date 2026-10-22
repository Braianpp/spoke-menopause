# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Q5: Recompute recurring genes for Subacute delirium → Dementia from top DWPC metapaths.

Endpoints
---------
- SOURCE: Disease:DOID:8645 (Subacute delirium)
- TARGET: Disease:DOID:1307 (Dementia)

Input metapaths (DWPC results)
-----------------------------
- ${SPOKE_DATA_ROOT}/biorouter/results/q5_dwpc_Subacute_delirium_Dementia.csv

Method
------
1) Filter metapaths to path_count > 0
2) Sort by dwpc desc; keep TOP_METAPATHS_LIMIT
3) For each metapath schema, enumerate up to MAX_PATHS paths constrained by edge attribute 'predicate'
4) Collect intermediate nodes with node_type == 'Gene'
5) Count recurrence across metapaths/paths

Output
------
- ${SPOKE_DATA_ROOT}/biorouter/results/q5_subacute_delirium_dementia_dwpc_recurring_genes.csv
"""
from __future__ import annotations
import os
import os
import pickle
from collections import Counter, defaultdict
from pathlib import Path
from typing import List, Tuple
import pandas as pd
SOURCE_NODE = 'Disease:DOID:8645'
TARGET_NODE = 'Disease:DOID:1307'
MAX_PATHS = int(os.environ.get('MAX_PATHS', '50'))
TOP_METAPATHS_LIMIT = int(os.environ.get('TOP_METAPATHS_LIMIT', '11'))
GFILE = os.environ.get('GFILE', '${SPOKE_DATA_ROOT}/V1_homosapiens_1N_9E.gpickle')
METAPATHS_CSV = os.environ.get('METAPATHS_CSV', '${SPOKE_DATA_ROOT}/biorouter/results/q5_dwpc_Subacute_delirium_Dementia.csv')
OUTCSV = os.environ.get('OUTCSV', '${SPOKE_DATA_ROOT}/biorouter/results/q5_subacute_delirium_dementia_dwpc_recurring_genes.csv')

def load_graph(gfile: str):
    p = Path(gfile)
    if not p.exists():
        raise SystemExit(f'Graph not found: {p}')
    with open(p, 'rb') as fh:
        return pickle.load(fh)

def _edge_dict_items(eds):
    """Yield (key, edata) pairs for both Graph and MultiGraph edge-data."""
    if not eds:
        return []
    if isinstance(eds, dict) and all((isinstance(v, dict) for v in eds.values())):
        return eds.items()
    return [(None, eds)]

def enumerate_paths_by_predicates(G, source: str, target: str, predicates: List[str], max_paths: int) -> List[List[str]]:
    """DFS enumeration constrained by a sequence of edge predicates.

    Returns a list of node-id paths (list[str]). Stops after max_paths successes.
    """
    paths: List[List[str]] = []
    stack = [source]

    def step(node: str, i: int):
        if len(paths) >= max_paths:
            return
        if i == len(predicates):
            if node == target:
                paths.append(stack.copy())
            return
        required = predicates[i]
        for nbr in G.neighbors(node):
            eds = G.get_edge_data(node, nbr)
            ok = False
            for _, edata in _edge_dict_items(eds):
                if edata.get('predicate') == required:
                    ok = True
                    break
            if not ok:
                continue
            stack.append(nbr)
            step(nbr, i + 1)
            stack.pop()
            if len(paths) >= max_paths:
                return
    step(source, 0)
    return paths

def main():
    mp_path = Path(METAPATHS_CSV)
    if not mp_path.exists():
        raise SystemExit(f'Missing metapaths CSV: {mp_path}')
    G = load_graph(GFILE)
    df_mp = pd.read_csv(mp_path)
    df_use = df_mp.copy()
    if 'path_count' in df_use.columns:
        df_use = df_use[df_use['path_count'] > 0]
    df_use = df_use.sort_values('dwpc', ascending=False).head(TOP_METAPATHS_LIMIT)
    gene_occ = Counter()
    gene_metapaths = defaultdict(set)
    for _, row in df_use.iterrows():
        schema = str(row['metapath']).split(' → ')
        paths = enumerate_paths_by_predicates(G, SOURCE_NODE, TARGET_NODE, schema, max_paths=MAX_PATHS)
        for p in paths:
            for n in p[1:-1]:
                nd = G.nodes[n]
                if nd.get('node_type') == 'Gene':
                    gene_occ[n] += 1
                    gene_metapaths[n].add(row['metapath'])
    if not gene_occ:
        raise SystemExit('No genes recovered (check metapaths, endpoints, or MAX_PATHS).')
    out = pd.DataFrame({'gene_id': list(gene_occ.keys()), 'total_occurrences': [gene_occ[g] for g in gene_occ.keys()], 'n_metapaths': [len(gene_metapaths[g]) for g in gene_occ.keys()]}).sort_values(['total_occurrences', 'n_metapaths', 'gene_id'], ascending=[False, False, True])
    out_path = Path(OUTCSV)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f'Saved: {out_path} (n_genes={len(out)})')
if __name__ == '__main__':
    main()
