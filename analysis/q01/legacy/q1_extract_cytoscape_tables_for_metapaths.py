# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Q1: Extract Cytoscape node/edge tables for specific metapaths (PM → Dementia).

Given a SPOKE graph and a list of metapath schemas (edge predicate sequences), this script:

1) Enumerates up to MAX_PATHS paths per metapath between:
   - SOURCE: Disease:DOID:10787 (Premature Menopause)
   - TARGET: Disease:DOID:1307  (Dementia)
2) Collects all nodes and edges that appear in any returned path.
3) Writes two CSVs suitable for Cytoscape import:
   - nodes table: node_id, node_type, name, n_paths
   - edges table: source, target, predicate, n_paths

Notes
-----
- The graph is treated as a (multi)graph; for each step we accept any parallel edge
  whose edata['predicate'] matches the required predicate.
- Edges are recorded in the *traversal direction* (source -> target) as the path is built.
- This is intended for *subgraph export* for visualization (not exhaustive enumeration).

Defaults:
- Graph: ${SPOKE_DATA_ROOT}/V1_homosapiens_1N_9E.gpickle
- Output dir: ${SPOKE_DATA_ROOT}/biorouter/results

Example
-------
python3 q1_extract_cytoscape_tables_for_metapaths.py   --max-paths 50   --out-prefix q1_pm_dementia_regulation
"""
from __future__ import annotations
import os
import argparse
import pickle
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
SOURCE_NODE = 'Disease:DOID:10787'
TARGET_NODE = 'Disease:DOID:1307'
DEFAULT_METAPATHS = ['ASSOCIATES_DaG → UPREGULATES_GPuG → ASSOCIATES_DaG', 'ASSOCIATES_DaG → UPREGULATES_KGuG → ASSOCIATES_DaG', 'ASSOCIATES_DaG → DOWNREGULATES_KGdG → ASSOCIATES_DaG', 'ASSOCIATES_DaG → DOWNREGULATES_GPdG → ASSOCIATES_DaG']

def load_graph(gfile: Path):
    if not gfile.exists():
        raise FileNotFoundError(f'Graph not found: {gfile}')
    with open(gfile, 'rb') as fh:
        return pickle.load(fh)

def enumerate_paths_by_predicates(G, source: str, target: str, predicates: List[str], max_paths: int) -> List[List[str]]:
    """DFS enumeration constrained by a sequence of edge predicates.

    Returns a list of node-id paths (each a list of node IDs).
    Stops after `max_paths` successful source->target paths.
    """
    paths: List[List[str]] = []
    stack: List[str] = [source]

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
            if not eds:
                continue
            if isinstance(eds, dict) and all((isinstance(v, dict) for v in eds.values())):
                items = eds.items()
            else:
                items = [(None, eds)]
            matched = False
            for _, edata in items:
                if edata.get('predicate') == required:
                    matched = True
                    break
            if not matched:
                continue
            stack.append(nbr)
            step(nbr, i + 1)
            stack.pop()
            if len(paths) >= max_paths:
                return
    step(source, 0)
    return paths

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--graph', default='${SPOKE_DATA_ROOT}/V1_homosapiens_1N_9E.gpickle', help='Path to SPOKE gpickle')
    ap.add_argument('--metapaths', nargs='*', default=DEFAULT_METAPATHS, help="Metapath strings using '→' separators")
    ap.add_argument('--max-paths', type=int, default=50, help='Max successful paths per metapath')
    ap.add_argument('--out-dir', default='${SPOKE_DATA_ROOT}/biorouter/results', help='Output directory')
    ap.add_argument('--out-prefix', default='q1_pm_dementia_selected_metapaths', help='Prefix for output CSV files')
    args = ap.parse_args()
    G = load_graph(Path(args.graph))
    metapath_schemas = [mp.split(' → ') for mp in args.metapaths]
    node_path_count = Counter()
    edge_path_count = Counter()
    total_paths = 0
    per_metapath_counts = []
    for mp_str, schema in zip(args.metapaths, metapath_schemas):
        paths = enumerate_paths_by_predicates(G, SOURCE_NODE, TARGET_NODE, schema, max_paths=args.max_paths)
        per_metapath_counts.append({'metapath': mp_str, 'n_paths': len(paths)})
        total_paths += len(paths)
        for p in paths:
            for n in set(p):
                node_path_count[n] += 1
            for i, pred in enumerate(schema):
                u = p[i]
                v = p[i + 1]
                edge_path_count[u, v, pred] += 1
    node_rows = []
    for node_id, n_paths in node_path_count.items():
        nd = G.nodes[node_id] if node_id in G else {}
        node_rows.append({'node_id': node_id, 'node_type': nd.get('node_type'), 'name': nd.get('name'), 'n_paths': int(n_paths)})
    nodes_df = pd.DataFrame(node_rows)
    type_order = {'Disease': 0, 'Gene': 1}
    nodes_df['_type_rank'] = nodes_df['node_type'].map(type_order).fillna(9)
    nodes_df = nodes_df.sort_values(['_type_rank', 'node_type', 'n_paths', 'node_id'], ascending=[True, True, False, True])
    nodes_df = nodes_df.drop(columns=['_type_rank'])
    edge_rows = []
    for (u, v, pred), n_paths in edge_path_count.items():
        edge_rows.append({'source': u, 'target': v, 'predicate': pred, 'n_paths': int(n_paths)})
    edges_df = pd.DataFrame(edge_rows).sort_values(['n_paths', 'predicate', 'source', 'target'], ascending=[False, True, True, True])
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    nodes_out = out_dir / f'{args.out_prefix}_nodes.csv'
    edges_out = out_dir / f'{args.out_prefix}_edges.csv'
    meta_out = out_dir / f'{args.out_prefix}_metapath_path_counts.csv'
    nodes_df.to_csv(nodes_out, index=False)
    edges_df.to_csv(edges_out, index=False)
    pd.DataFrame(per_metapath_counts).to_csv(meta_out, index=False)
    print(f'Saved nodes: {nodes_out} (n={len(nodes_df)})')
    print(f'Saved edges: {edges_out} (n={len(edges_df)})')
    print(f'Saved metapath counts: {meta_out}')
    print(f'Total enumerated paths (capped): {total_paths}')
if __name__ == '__main__':
    main()
