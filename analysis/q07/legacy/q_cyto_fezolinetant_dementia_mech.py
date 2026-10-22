# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Fezolinetant -> Dementia, strict mechanistic (option C) for Cytoscape.

Layer 1 (direct mechanistic): Compound -BINDS- Protein <-ENCODES- Gene -ASSOCIATES_DaG- Disease(dementia)
Layer 2 (target-mediated repurposing via NEUROKININ receptors ONLY):
   Compound -BINDS- Protein(NK1R/NK3R) -INTERACTS_PiC|BINDS_CbP- Compound2 -TREATS_CtD|IN_CLINICAL_TRIALS_FOR_CictD- Disease(dementia)
   (CYP/KCNH2-mediated bridges are excluded -> pharmacokinetic, not mechanistic)
Dedup: if a node pair has multiple edges and one is TREATS_CtD, keep only TREATS_CtD.
"""
import os
import sys, csv, os
sys.path.insert(0, '${SPOKE_DATA_ROOT}/biorouter/scripts')
from spoke_utils import open_spoke
OUTDIR = '${SPOKE_DATA_ROOT}/biorouter/exports'
CID = 'inchikey:PPSNFPASKFYPMN-SECBINFHSA-N'
TERMS = ['dementia', 'alzheimer', 'neurodegener', 'cognit', 'frontotemporal', 'lewy', 'tauopath', 'mild cognitive']
NK_GENES = ['TACR1', 'TACR3', 'TACR2']
nodes = {}
edges = {}

def AN(i, l, t):
    if i not in nodes:
        nodes[i] = (l, t)

def AE(s, t, inter, et):
    edges[s, t, inter] = et
drv = open_spoke()
with drv.session() as s:
    W = 'any(x IN $T WHERE toLower(d.name) CONTAINS x)'
    q1 = 'MATCH (c:Compound {identifier:$id})-[:BINDS_CbP]-(p:Protein)<-[:ENCODES_GeP]-(g:Gene)-[:ASSOCIATES_DaG]-(d:Disease) WHERE ' + W + ' RETURN c.identifier AS ci,c.name AS cn,p.identifier AS pi,p.name AS pn,g.name AS gn,d.identifier AS di,d.name AS dn'
    for r in s.run(q1, id=CID, T=TERMS).data():
        AN(r['ci'], r['cn'], 'Compound (query)')
        AN(r['pi'], r['pn'], 'Protein (shared target)')
        AN('gene:' + r['gn'], r['gn'], 'Gene')
        AN(r['di'], r['dn'], 'Disease (dementia)')
        AE(r['ci'], r['pi'], 'BINDS_CbP', 'BINDS_CbP')
        AE('gene:' + r['gn'], r['pi'], 'ENCODES_GeP', 'ENCODES_GeP')
        AE('gene:' + r['gn'], r['di'], 'ASSOCIATES_DaG', 'ASSOCIATES_DaG')
    q2 = 'MATCH (c:Compound {identifier:$id})-[:BINDS_CbP]-(p:Protein)<-[:ENCODES_GeP]-(gp:Gene) WHERE gp.name IN $NK MATCH (p)-[e2]-(c2:Compound)-[e3]-(d:Disease) WHERE ' + W + " AND type(e2) IN ['INTERACTS_PiC','BINDS_CbP'] AND type(e3) IN ['IN_CLINICAL_TRIALS_FOR_CictD','TREATS_CtD'] AND c2.identifier<>$id RETURN c.identifier AS ci,c.name AS cn,p.identifier AS pi,p.name AS pn,type(e2) AS e2,c2.identifier AS c2i,c2.name AS c2n,type(e3) AS e3,d.identifier AS di,d.name AS dn"
    for r in s.run(q2, id=CID, T=TERMS, NK=NK_GENES).data():
        AN(r['ci'], r['cn'], 'Compound (query)')
        AN(r['pi'], r['pn'], 'Protein (shared target)')
        AN(r['c2i'], r['c2n'], 'Compound (bridge)')
        AN(r['di'], r['dn'], 'Disease (dementia)')
        AE(r['ci'], r['pi'], 'BINDS_CbP', 'BINDS_CbP')
        AE(r['pi'], r['c2i'], r['e2'], r['e2'])
        AE(r['c2i'], r['di'], r['e3'], r['e3'])
drv.close()
from collections import defaultdict
pair_edges = defaultdict(list)
for (s, t, inter), et in edges.items():
    pair_edges[frozenset((s, t))].append((s, t, inter, et))
final_edges = {}
for pair, lst in pair_edges.items():
    treats = [e for e in lst if e[3] == 'TREATS_CtD']
    keep = treats if treats else lst
    for s, t, inter, et in keep:
        final_edges[s, t, inter] = et
removed = len(edges) - len(final_edges)
print('EDGES_BEFORE_DEDUP:', len(edges), ' AFTER:', len(final_edges), ' REMOVED:', removed)
keep_ids = set(nodes)
disease = set((i for i, (l, t) in nodes.items() if t == 'Disease (dementia)'))
bridge = set((i for i, (l, t) in nodes.items() if t == 'Compound (bridge)'))
reach = set()
for s, t, inter in final_edges:
    if s in bridge and t in disease:
        reach.add(s)
    if t in bridge and s in disease:
        reach.add(t)
orphans = bridge - reach
nodes = {i: v for i, v in nodes.items() if i not in orphans}
final_edges = {k: v for k, v in final_edges.items() if k[0] in nodes and k[1] in nodes}
os.makedirs(OUTDIR, exist_ok=True)
nf = OUTDIR + '/fezo_dementia_mech_nodes.csv'
ef = OUTDIR + '/fezo_dementia_mech_edges.csv'
with open(nf, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['id', 'label', 'node_type'])
    for i, (l, t) in nodes.items():
        w.writerow([i, l, t])
with open(ef, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['source', 'target', 'interaction', 'edge_type', 'directed'])
    for (s, t, inter), et in final_edges.items():
        w.writerow([s, t, inter, et, 'true'])
from collections import Counter
print('NODES:', len(nodes), dict(Counter((t for _, t in nodes.values()))))
print('EDGES:', len(final_edges), dict(Counter(final_edges.values())))
print('PROTEINS:', sorted((l for l, t in nodes.values() if t == 'Protein (shared target)')))
print('BRIDGES:', sorted((l for l, t in nodes.values() if t == 'Compound (bridge)')))
print('DISEASES:', sorted((l for l, t in nodes.values() if t == 'Disease (dementia)')))
print('GENES:', sorted((l for l, t in nodes.values() if t == 'Gene')))
print(nf)
print(ef)
