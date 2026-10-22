# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Fezolinetant -> Breast cancer pathway (WP4262), fully connected metapath, for Cytoscape.

Metapath (from original figure):
  Fezolinetant -[BINDS_CbP]- P0(target) -[INTERACTS_PiP]- P1(pathway protein)
    <-[ENCODES_GeP]- Gene -[PARTICIPATES_GpPW]- Pathway(WP4262 Breast cancer pathway)
Only keeps P1/Gene that actually participate in WP4262, so the pathway node is connected.
"""
import os
import sys, csv, os
sys.path.insert(0, '${SPOKE_DATA_ROOT}/biorouter/scripts')
from spoke_utils import open_spoke
OUTDIR = '${SPOKE_DATA_ROOT}/biorouter/exports'
CID = 'inchikey:PPSNFPASKFYPMN-SECBINFHSA-N'
PW = 'WP4262'
nodes = {}
edges = {}

def AN(i, l, t):
    if i not in nodes:
        nodes[i] = (l, t)

def AE(s, t, inter):
    edges[s, t, inter] = inter
drv = open_spoke()
with drv.session() as s:
    q = 'MATCH (c:Compound {identifier:$id})-[:BINDS_CbP]-(p0:Protein) MATCH (p0)-[:INTERACTS_PiP]-(p1:Protein)<-[:ENCODES_GeP]-(g:Gene)-[:PARTICIPATES_GpPW]-(pw:Pathway {identifier:$pw}) RETURN DISTINCT c.identifier AS ci,c.name AS cn, p0.identifier AS p0i,p0.name AS p0n, p1.identifier AS p1i,p1.name AS p1n, g.identifier AS gi,g.name AS gn, pw.identifier AS pwi,pw.name AS pwn'
    rows = s.run(q, id=CID, pw=PW).data()
    print('METAPATH_ROWS:', len(rows))
    for r in rows:
        AN('Compound:' + r['ci'], r['cn'], 'Compound')
        AN('Protein:' + r['p0i'], r['p0n'], 'Protein')
        AN('Protein:' + r['p1i'], r['p1n'], 'Protein')
        AN('Gene:' + str(r['gi']), r['gn'], 'Gene')
        AN('Pathway:' + r['pwi'], r['pwn'], 'Pathway')
        AE('Compound:' + r['ci'], 'Protein:' + r['p0i'], 'BINDS_CbP')
        AE('Protein:' + r['p0i'], 'Protein:' + r['p1i'], 'INTERACTS_PiP')
        AE('Gene:' + str(r['gi']), 'Protein:' + r['p1i'], 'ENCODES_GeP')
        AE('Gene:' + str(r['gi']), 'Pathway:' + r['pwi'], 'PARTICIPATES_GpPW')
drv.close()
os.makedirs(OUTDIR, exist_ok=True)
nf = OUTDIR + '/fezo_breast_connected_nodes.csv'
ef = OUTDIR + '/fezo_breast_connected_edges.csv'
with open(nf, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['id', 'label', 'node_type'])
    for i, (l, t) in nodes.items():
        w.writerow([i, l, t])
with open(ef, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['source', 'target', 'interaction', 'edge_type', 'directed'])
    for s, t, inter in edges:
        w.writerow([s, t, inter, inter, 'true'])
from collections import Counter
print('NODES:', len(nodes), dict(Counter((t for _, t in nodes.values()))))
print('EDGES:', len(edges), dict(Counter(edges.values())))
print('TARGETS(P0):', sorted(set((l for i, (l, t) in nodes.items() if t == 'Protein' and any((i == 'Protein:' + r['p0i'] for r in rows))))))
print('GENES:', sorted((l for l, t in nodes.values() if t == 'Gene')))
print('PATHWAY:', [l for l, t in nodes.values() if t == 'Pathway'])
print(nf)
print(ef)
