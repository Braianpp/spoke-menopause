# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
import os, sys, pickle, time
from pathlib import Path
from typing import List, Dict, Set, Tuple, Optional
from collections import defaultdict
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import sparse
try:
    import networkx as nx
except Exception:
    nx = None
BASE = Path('${SPOKE_DATA_ROOT}')
OUTDIR = BASE / 'biorouter' / 'results'
OUTDIR.mkdir(parents=True, exist_ok=True)
LOGDIR = BASE / 'biorouter' / 'logs'
LOGDIR.mkdir(parents=True, exist_ok=True)
CANDS = [BASE / 'spoke_graph.gpickle', Path('${SPOKE_DATA_ROOT}/V1_homosapiens_1N_9E.gpickle'), Path('${SPOKE_DATA_ROOT}/extract_spoke/spokeML/extract_spoke/V1_homosapiens_1N_9E.gpickle')]
G = None
for gfile in CANDS:
    try:
        if nx is not None and hasattr(nx, 'read_gpickle'):
            G = nx.read_gpickle(gfile)
        else:
            raise AttributeError('no read_gpickle')
    except Exception:
        G = None
        if gfile.exists():
            try:
                with open(gfile, 'rb') as fh:
                    G = pickle.load(fh)
            except Exception:
                G = None
    if G is not None:
        print(f'[q5_dwpc_compute_pairs] Loaded graph from {gfile}')
        break
if G is None:
    print('[q5_dwpc_compute_pairs] ERROR: Could not load G from known paths', file=sys.stderr)
    sys.exit(2)

class NetworkDWPC:

    def __init__(self, G: nx.MultiGraph, damping_factor: float=-0.4):
        self.G = G
        self.damping_factor = damping_factor
        self.metapath_finder = MetapathFinder(G)
        self.node_to_idx = {node: idx for idx, node in enumerate(G.nodes())}
        self.idx_to_node = {idx: node for node, idx in self.node_to_idx.items()}
        self.edge_matrices: Dict[str, sparse.csr_matrix] = {}

    def find_metapaths(self, source_type: str, target_type: str, max_length: int=2) -> List[List[str]]:
        return self.metapath_finder.find_metapaths(source_type, target_type, max_length)

    def _get_edge_matrix(self, et: str) -> Optional[sparse.csr_matrix]:
        if et in self.edge_matrices:
            return self.edge_matrices[et]
        n = len(self.G)
        M = sparse.lil_matrix((n, n))
        found = False
        for u, v, data in self.G.edges(data=True):
            if data.get('predicate') == et:
                ui = self.node_to_idx[u]
                vi = self.node_to_idx[v]
                M[ui, vi] = 1
                if not self.G.is_directed():
                    M[vi, ui] = 1
                found = True
        if not found:
            return None
        C = M.tocsr()
        self.edge_matrices[et] = C
        return C

    def _create_degree_matrices(self, metapath: List[str]) -> List[Tuple[sparse.dia_matrix, sparse.dia_matrix]]:
        degree_mats = []
        for et in metapath:
            A = self.edge_matrices.get(et)
            if A is None:
                out_deg = in_deg = np.ones(len(self.G))
            else:
                out_deg = np.array(A.sum(axis=1)).ravel()
                in_deg = np.array(A.sum(axis=0)).ravel()
                out_deg[out_deg == 0] = 1
                in_deg[in_deg == 0] = 1
            degree_mats.append((sparse.diags(np.power(out_deg, self.damping_factor)), sparse.diags(np.power(in_deg, self.damping_factor))))
        return degree_mats

    def compute_dwpc(self, source: str, target: str, metapath: List[str]) -> float:
        if source not in self.node_to_idx or target not in self.node_to_idx:
            return 0.0
        mats = []
        for et in metapath:
            A = self._get_edge_matrix(et)
            if A is None:
                return 0.0
            mats.append(A)
        s = self.node_to_idx[source]
        t = self.node_to_idx[target]
        deg = self._create_degree_matrices(metapath)
        result = mats[0]
        result = deg[0][0].dot(result)
        for i in range(len(metapath) - 1):
            result = result.dot(deg[i][1])
            result = result.dot(deg[i + 1][0])
            result = result.dot(mats[i + 1])
        result = result.dot(deg[-1][1])
        return float(result[s, t])

    def compute_simple_path_count(self, source: str, target: str, metapath: List[str]) -> int:
        if source not in self.node_to_idx or target not in self.node_to_idx:
            return 0
        mats = []
        for et in metapath:
            A = self._get_edge_matrix(et)
            if A is None:
                return 0
            mats.append(A)
        s = self.node_to_idx[source]
        t = self.node_to_idx[target]
        result = mats[0]
        for A in mats[1:]:
            result = result.dot(A)
        return int(result[s, t])

class MetapathFinder:

    def __init__(self, G: nx.MultiGraph):
        self.G = G
        self._cache_network_structure()

    def _cache_network_structure(self):
        self.node_types = set(nx.get_node_attributes(self.G, 'node_type').values())
        self.valid_transitions = defaultdict(set)
        self.edge_types = set()
        for u, v, data in self.G.edges(data=True):
            et = data.get('predicate')
            if et is None:
                continue
            ut = self.G.nodes[u].get('node_type')
            vt = self.G.nodes[v].get('node_type')
            if ut is None or vt is None:
                continue
            self.valid_transitions[ut, vt].add(et)
            self.edge_types.add(et)
            if not self.G.is_directed():
                self.valid_transitions[vt, ut].add(et)

    def _extend(self, cur_path: List[str], cur_type: str, tgt_type: str, length_left: int, out: List[List[str]]):
        if length_left == 0:
            if cur_type == tgt_type:
                out.append(list(cur_path))
            return
        for next_type in self.node_types:
            ets = self.valid_transitions.get((cur_type, next_type), set())
            for et in ets:
                cur_path.append(et)
                self._extend(cur_path, next_type, tgt_type, length_left - 1, out)
                cur_path.pop()

    def find_metapaths(self, source_type: str, target_type: str, max_length: int=2) -> List[List[str]]:
        mets: List[List[str]] = []
        for L in range(1, max_length + 1):
            self._extend([], source_type, target_type, L, mets)
        return mets
pairs = [('Disease:DOID:10787', 'Premature menopause', 'Disease:DOID:1307', 'Dementia'), ('Disease:DOID:8645', 'Subacute delirium', 'Disease:DOID:1307', 'Dementia')]
calc = NetworkDWPC(G)
MAX_LEN = int(os.environ.get('Q5_MAX_LENGTH', '3'))
all_rows = []
for sid, sname, tid, tname in pairs:
    if sid not in G.nodes or tid not in G.nodes:
        print(f'[q5_dwpc_compute_pairs] WARNING: missing node for {sname}→{tname}')
        continue
    st = G.nodes[sid].get('node_type')
    tt = G.nodes[tid].get('node_type')
    mets = calc.find_metapaths(st, tt, max_length=MAX_LEN)
    print(f'[q5_dwpc_compute_pairs] {sname} → {tname}: {len(mets)} metapaths (<= {MAX_LEN})')
    for mp in tqdm(mets, desc=f'{sname}→{tname}'):
        dwpc = calc.compute_dwpc(sid, tid, mp)
        pcnt = calc.compute_simple_path_count(sid, tid, mp)
        all_rows.append({'pair': f'{sname} → {tname}', 'source_id': sid, 'source_name': sname, 'target_id': tid, 'target_name': tname, 'metapath': ' → '.join(mp), 'length': len(mp), 'dwpc': dwpc, 'path_count': pcnt})
res = pd.DataFrame(all_rows)
if res.empty:
    print('[q5_dwpc_compute_pairs] WARNING: no results computed')
res.to_csv(OUTDIR / 'q5_dwpc_pairs.csv', index=False)
for pair_label, dfp in res.groupby('pair'):
    safe = pair_label.replace(' → ', '_').replace(' ', '_')
    dfp.sort_values('dwpc', ascending=False).to_csv(OUTDIR / f'q5_dwpc_{safe}.csv', index=False)
for pair_label, dfp in res.groupby('pair'):
    nz = dfp[dfp['dwpc'] > 0].copy()
    print('\n' + '=' * 60)
    print(pair_label)
    print(f"Active metapaths: {len(nz)}; cumulative DWPC = {nz['dwpc'].sum():.6e}")
    print(nz.sort_values('dwpc', ascending=False).head(10)[['metapath', 'dwpc', 'path_count']].to_string(index=False))
print(f"\nSaved: {OUTDIR / 'q5_dwpc_pairs.csv'}")
