# Portable replay adapted from matched_null.py; original preserved alongside.
# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Reproduce the historical real extraction and apply it identically to 1500 pairs.
Exact ordered DFS for schemas <=3 edges, with reverse reachability pruning.
No expansion cap, no randomized neighbor order, no swallowed failures.
"""
import os
from pathlib import Path
from collections import Counter, defaultdict
import pickle, random, json, hashlib, time, platform
import pandas as pd
import numpy as np
import networkx as nx
import argparse
from spoke_menopause.paper import load_tables, records
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--graph',type=Path,required=True)
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args()
ROOT=Path(__file__).resolve().parents[2]
OUT=args.out
OUT.mkdir(parents=True,exist_ok=True)
GF=args.graph
MF=ROOT/'data/q06/selected_schemas.csv'
S, T = ('Disease:DOID:8645', 'Disease:DOID:1307')
N = 1500
CAP = 50
SEED = 42

def status(stage, **kw):
    obj = {'stage': stage, 'utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), **kw}
    (OUT / 'status.json').write_text(json.dumps(obj, indent=2))
    print(json.dumps(obj), flush=True)
status('loading graph')
with GF.open('rb') as f:
    G = pickle.load(f)
assert not G.is_directed()
mp = pd.read_csv(MF)
mp = mp[mp.path_count > 0].sort_values('dwpc', ascending=False).head(11)
schemas = [x.split(' → ') for x in mp.metapath]
assert all((1 <= len(x) <= 3 for x in schemas))
mp.to_csv(OUT / 'selected_schemas.csv', index=False)
cache = {}
pred_cache = {}

def adj(n, p):
    key = (n, p)
    if key not in cache:
        ns = []
        for nb, eds in G[n].items():
            vals = eds.values() if G.is_multigraph() else [eds]
            if any((e.get('predicate') == p for e in vals)):
                ns.append(nb)
        cache[key] = (ns, set(ns), {x: i for i, x in enumerate(ns)})
    return cache[key]

def intersect_ordered(n, p, allowed):
    ns, st, order = adj(n, p)
    if len(ns) <= len(allowed):
        return [x for x in ns if x in allowed]
    return sorted((x for x in allowed if x in st), key=order.__getitem__)

def paths(s, t, ps):
    if len(ps) == 1:
        return [[s, t]] if t in adj(s, ps[0])[1] else []
    last = adj(t, ps[-1])[1]
    if len(ps) == 2:
        return [[s, a, t] for a in intersect_ordered(s, ps[0], last)[:CAP]]
    out = []
    for a in adj(s, ps[0])[0]:
        for b in intersect_ordered(a, ps[1], last):
            out.append([s, a, b, t])
            if len(out) >= CAP:
                return out
    return out

def extract(s, t):
    ctr = Counter()
    rec = defaultdict(set)
    allpaths = []
    for idx, ps in enumerate(schemas):
        found = paths(s, t, ps)
        for path in found:
            allpaths.append({'schema_index': idx, 'nodes': path})
            for n in path[1:-1]:
                if G.nodes[n].get('node_type') == 'Gene':
                    ctr[n] += 1
                    rec[n].add(idx)
    return (ctr, rec, allpaths)
real, rec, rpaths = extract(S, T)
old = pd.DataFrame(records(load_tables(ROOT),25)).rename(columns={'Gene identifier':'gene_id','Observed recurrence':'total_occurrences','Metapaths containing gene':'n_metapaths'}).set_index('gene_id')
assert real == Counter(old.total_occurrences.to_dict()), 'Real counts changed: stop and investigate'
assert all((len(rec[g]) == int(old.loc[g, 'n_metapaths']) for g in real))
(OUT / 'real_paths.json').write_text(json.dumps(rpaths, indent=2))
status('real reproduced', genes=len(real), paths=len(rpaths))

def preds(n):
    if n not in pred_cache:
        out = set()
        for es in G[n].values():
            for e in es.values() if G.is_multigraph() else [es]:
                if e.get('predicate'):
                    out.add(e['predicate'])
        pred_cache[n] = out
    return pred_cache[n]
diseases = [n for n, d in G.nodes(data=True) if str(d.get('node_type', '')).lower() == 'disease']
srcp = {p[0] for p in schemas}
tgtp = {p[-1] for p in schemas}
rng = random.Random(SEED)

def draw(req):
    for _ in range(100000):
        n = rng.choice(diseases)
        if n not in (S, T) and preds(n) & req:
            return n
    raise RuntimeError('No compatible disease sampled')
pairs = list(pd.read_csv(ROOT/'data/q06/sampled_pairs.csv').itertuples(index=False,name=None))
assert len(pairs)==N
pd.DataFrame(pairs, columns=['source', 'target']).to_csv(OUT / 'sampled_pairs.csv', index=False)
null = []
yields = []
with (OUT / 'null_gene_counts.jsonl').open('w') as f:
    for i, (a, b) in enumerate(pairs, 1):
        c, _, pp = extract(a, b)
        null.append(c)
        yields.append(len(c))
        f.write(json.dumps({'pair': i, 'source': a, 'target': b, 'counts': c, 'unique_genes': len(c), 'paths': len(pp)}) + '\n')
        if i % 50 == 0:
            f.flush()
            status('null running', completed=i, total=N)
rows = []
for g, c in real.items():
    vals = np.array([d.get(g, 0) for d in null])
    nd = G.nodes[g]
    rows.append({'gene_id': g, 'gene_symbol': nd.get('name', nd.get('symbol', g)), 'real_count': c, 'n_metapaths': len(rec[g]), 'bg_mean': float(vals.mean()), 'bg_std': float(vals.std()), 'bg_max': int(vals.max()), 'n_geq': int((vals >= c).sum()), 'empirical_pval': float((1 + (vals >= c).sum()) / (N + 1))})
df = pd.DataFrame(rows).sort_values('empirical_pval')
p = df.empirical_pval.to_numpy()
q = np.minimum.accumulate((p * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
df['adjusted_pval'] = np.minimum(q, 1)
df['significant_bh'] = df.adjusted_pval < 0.05
df.to_csv(OUT / 'matched_gene_null.csv', index=False)
pd.DataFrame({'pair': range(1, N + 1), 'bg_unique_genes': yields}).to_csv(OUT / 'matched_null_yield.csv', index=False)
info = {'n_pairs': N, 'n_unique_pairs': len(set(pairs)), 'seed': SEED, 'path_cap_per_schema': CAP, 'expansion_cap': None, 'schemas': len(schemas), 'real_genes': len(real), 'real_counts_match_archived': True, 'null_gene_yield_mean': float(np.mean(yields)), 'null_gene_yield_max': max(yields), 'yield_n_geq': sum((x >= len(real) for x in yields)), 'yield_empirical_p': (1 + sum((x >= len(real) for x in yields))) / (N + 1), 'significant_genes': int(df.significant_bh.sum()), 'max_adjusted_pval': float(df.adjusted_pval.max()), 'graph_nodes': len(G), 'graph_edges': G.number_of_edges(), 'python': platform.python_version(), 'networkx': nx.__version__, 'graph_filename': GF.name, 'graph_mtime': GF.stat().st_mtime, 'schema_sha256': hashlib.sha256(MF.read_bytes()).hexdigest(), 'method': 'Same ordered enumeration for real and null; reverse pruning preserves original DFS order; replay of saved predicate-compatible ordered disease pairs originally sampled with replacement; fixed observed schemas; not degree matched.'}
(OUT / 'summary.json').write_text(json.dumps(info, indent=2))
status('complete', **info)
