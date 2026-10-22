# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""
Elinzanetant -> Dementia network (<=3 hops) for Cytoscape.

Path pattern (data-driven from shortestPath diagnostics):
  Compound(Elinzanetant) -[BINDS_CbP]- Protein -[INTERACTS_PiC|BINDS_CbP]- Compound2
    -[IN_CLINICAL_TRIALS_FOR_CictD|TREATS_CtD]- Disease(dementia-family)

This is a shared-target / drug-repurposing bridge (not a gene->disease mechanism):
Elinzanetant shares NK1R/NK3R targets with compounds trialed for / treating dementia.

Outputs (Cytoscape import format):
  elinza_dementia_3hop_nodes.csv  : id,label,node_type
  elinza_dementia_3hop_edges.csv  : source,target,interaction,edge_type,directed
"""
import os
import sys, csv, os
sys.path.insert(0, '${SPOKE_DATA_ROOT}/biorouter/scripts')
from spoke_utils import open_spoke
OUTDIR = '${SPOKE_DATA_ROOT}/biorouter/exports'
os.makedirs(OUTDIR, exist_ok=True)
CID = 'inchikey:DWRIJNIPBUFCQS-DQEYMECFSA-N'
TERMS = ['dementia', 'alzheimer', 'neurodegener', 'cognit', 'frontotemporal', 'lewy', 'tauopath', 'mild cognitive']
CY = "\nMATCH (c:Compound {identifier:$id})-[e1:BINDS_CbP]-(p:Protein)-[e2]-(c2:Compound)-[e3]-(d:Disease)\nWHERE any(t IN $TERMS WHERE toLower(d.name) CONTAINS t)\n  AND type(e2) IN ['INTERACTS_PiC','BINDS_CbP']\n  AND type(e3) IN ['IN_CLINICAL_TRIALS_FOR_CictD','TREATS_CtD']\n  AND c2.identifier <> $id\nRETURN c.identifier AS c_id, c.name AS c_name,\n       p.identifier AS p_id, p.name AS p_name,\n       type(e2) AS e2, c2.identifier AS c2_id, c2.name AS c2_name,\n       type(e3) AS e3, d.identifier AS d_id, d.name AS d_name\n"
nodes = {}
edges = {}

def add_node(nid, label, ntype):
    if nid not in nodes:
        nodes[nid] = (label, ntype)

def add_edge(src, tgt, interaction, etype):
    edges[src, tgt, interaction] = etype
drv = open_spoke()
with drv.session() as s:
    rows = s.run(CY, id=CID, TERMS=TERMS).data()
    for r in rows:
        add_node(r['c_id'], r['c_name'], 'Compound (query)')
        add_node(r['p_id'], r['p_name'], 'Protein (shared target)')
        add_node(r['c2_id'], r['c2_name'], 'Compound (bridge)')
        add_node(r['d_id'], r['d_name'], 'Disease (dementia)')
        add_edge(r['c_id'], r['p_id'], 'BINDS_CbP', 'BINDS_CbP')
        add_edge(r['p_id'], r['c2_id'], r['e2'], r['e2'])
        add_edge(r['c2_id'], r['d_id'], r['e3'], r['e3'])
drv.close()
nf = os.path.join(OUTDIR, 'elinza_dementia_3hop_nodes.csv')
ef = os.path.join(OUTDIR, 'elinza_dementia_3hop_edges.csv')
with open(nf, 'w', newline='') as f:
    wc = csv.writer(f)
    wc.writerow(['id', 'label', 'node_type'])
    for nid, (label, ntype) in nodes.items():
        wc.writerow([nid, label, ntype])
with open(ef, 'w', newline='') as f:
    wc = csv.writer(f)
    wc.writerow(['source', 'target', 'interaction', 'edge_type', 'directed'])
    for (src, tgt, inter), etype in edges.items():
        wc.writerow([src, tgt, inter, etype, 'true'])
print('NODES:', len(nodes), '-> ', nf)
print('EDGES:', len(edges), '-> ', ef)
from collections import Counter
print('NODE_TYPES:', dict(Counter((t for _, t in nodes.values()))))
print('EDGE_TYPES:', dict(Counter(edges.values())))
print('DISEASES:', [l for l, t in nodes.values() if t == 'Disease (dementia)'])
print('BRIDGE_COMPOUNDS:', [l for l, t in nodes.values() if t == 'Compound (bridge)'])
