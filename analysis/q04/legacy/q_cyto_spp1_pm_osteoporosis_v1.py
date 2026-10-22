# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""
q_cyto_spp1_pm_osteoporosis_v1.py

Build a Cytoscape-style network figure showing SPP1 (gene; protein osteopontin,
UniProt P10451) as the candidate MEDIATOR linking premature menopause (POI) and
osteoporosis, with the bone-brain (dementia) extension. For manuscript Q4 Fig. 5b.

Method (matches q_cyto_pm_dementia_bridge_genes_v1.py conventions):
- Query SPOKE dev graph on leo via the neo4j driver + spoke_utils credentials.
- Extract only edges that lie on a premature-menopause <-> SPP1 <-> osteoporosis
  path, plus the SPP1 <-> dementia bone-brain extension.

Metapaths used (all edges are REAL SPOKE query results, nothing fabricated):
  A. SPP1 -ASSOCIATES_DaG-> osteoporosis                      (direct)
  B. premature menopause -DaG- g2 -GpPW- pathway -GpPW- SPP1  (pathway bridge to PM)
  C. SPP1 protein (P10451) -INTERACTS_PiP- partner protein
       <-ENCODES_GeP- gene -ASSOCIATES_DaG-> dementia/Alzheimer (bone-brain PPI)
  D. dementia/Alzheimer -DaG- g -GpPW- shared bone pathway -GpPW- SPP1
       restricted to bone/inflammatory pathways (osteoporosis-dementia bridge)

Outputs (written to leo exports, then scp'd back locally by the caller):
  <OUT_DIR>/spp1_pm_osteoporosis_nodes.csv
  <OUT_DIR>/spp1_pm_osteoporosis_edges.csv
"""
import os
import os, sys
import pandas as pd
from neo4j import GraphDatabase
sys.path.append('${SPOKE_DATA_ROOT}/biorouter/scripts')
sys.path.append('${SPOKE_DATA_ROOT}')
from spoke_utils import SPOKE_DEV_URI, SPOKE_User, SPOKE_Password
OUT_DIR = os.environ.get('CYTO_OUT', '${SPOKE_DATA_ROOT}/biorouter/exports')
os.makedirs(OUT_DIR, exist_ok=True)
PM_ID = 'DOID:10787'
OSTEO_ID = 'DOID:11476'
SPP1 = 'spp1'
OPN = 'P10451'
KEEP_DEM_PARTNERS = {'IL6', 'IL1B', 'TNF', 'TREM2', 'APOE'}
driver = GraphDatabase.driver(SPOKE_DEV_URI, auth=(SPOKE_User, SPOKE_Password))
nodes = {}
edges = {}

def add_node(ntype, ident, name, role=''):
    nid = f'{ntype}:{ident}'
    if nid not in nodes:
        nodes[nid] = {'id': nid, 'name': name, 'type': ntype, 'role': role}
    elif role and (not nodes[nid]['role']):
        nodes[nid]['role'] = role
    return nid

def add_edge(src, tgt, interaction, abbr):
    key = (src, tgt, interaction)
    if key not in edges:
        edges[key] = {'source': src, 'target': tgt, 'interaction': interaction, 'abbr': abbr}
with driver.session() as S:
    spp1_row = S.run('MATCH (g:Gene) WHERE toLower(g.name)=$n RETURN g.identifier AS id, g.name AS name', n=SPP1).single()
    SPP1_ID = spp1_row['id']
    add_node('Gene', SPP1_ID, spp1_row['name'], role='mediator')
    add_node('Protein', OPN, 'Osteopontin (SPP1)', role='mediator')
    add_edge(f'Gene:{SPP1_ID}', f'Protein:{OPN}', 'ENCODES_GeP', 'GeP')
    pm_row = S.run('MATCH (d:Disease {identifier:$i}) RETURN d.name AS name', i=PM_ID).single()
    add_node('Disease', PM_ID, pm_row['name'], role='Premature Menopause')
    osteo_row = S.run('MATCH (d:Disease {identifier:$i}) RETURN d.name AS name', i=OSTEO_ID).single()
    add_node('Disease', OSTEO_ID, osteo_row['name'], role='Osteoporosis')
    for r in S.run('\n        MATCH (g:Gene {identifier:$g})-[:ASSOCIATES_DaG]-(d:Disease {identifier:$o})\n        RETURN d.identifier AS id\n    ', g=SPP1_ID, o=OSTEO_ID).data():
        add_edge(f'Gene:{SPP1_ID}', f'Disease:{OSTEO_ID}', 'ASSOCIATES_DaG', 'DaG')
    for r in S.run('\n        MATCH (pm:Disease {identifier:$pm})-[:ASSOCIATES_DaG]-(g2:Gene)\n              -[:PARTICIPATES_GpPW]-(pw:Pathway)-[:PARTICIPATES_GpPW]-(spp1:Gene {identifier:$g})\n        RETURN DISTINCT g2.identifier AS g2_id, g2.name AS g2_name,\n                        pw.identifier AS pw_id, pw.name AS pw_name\n    ', pm=PM_ID, g=SPP1_ID).data():
        g2 = add_node('Gene', r['g2_id'], r['g2_name'], role='PM_gene')
        pw = add_node('Pathway', r['pw_id'], r['pw_name'])
        add_edge(f'Disease:{PM_ID}', g2, 'ASSOCIATES_DaG', 'DaG')
        add_edge(g2, pw, 'PARTICIPATES_GpPW', 'GpPW')
        add_edge(pw, f'Gene:{SPP1_ID}', 'PARTICIPATES_GpPW', 'GpPW')
    for r in S.run("\n        MATCH (op:Protein {identifier:$op})-[:INTERACTS_PiP]-(p2:Protein)\n              -[:ENCODES_GeP]-(g2:Gene)-[:ASSOCIATES_DaG]-(dem:Disease)\n        WHERE (toLower(dem.name) CONTAINS 'dementia' OR toLower(dem.name) CONTAINS 'alzheimer')\n          AND g2.name IN $keep\n        RETURN DISTINCT dem.identifier AS dem_id, dem.name AS dem_name,\n                        g2.identifier AS g2_id, g2.name AS g2_name,\n                        p2.identifier AS p2_id, p2.name AS p2_name\n    ", op=OPN, keep=list(KEEP_DEM_PARTNERS)).data():
        dem = add_node('Disease', r['dem_id'], r['dem_name'], role='Dementia')
        g2 = add_node('Gene', r['g2_id'], r['g2_name'], role='inflammatory')
        p2 = add_node('Protein', r['p2_id'], r['p2_name'])
        add_edge(f'Protein:{OPN}', p2, 'INTERACTS_PiP', 'PiP')
        add_edge(p2, g2, 'ENCODES_GeP', 'GeP')
        add_edge(g2, dem, 'ASSOCIATES_DaG', 'DaG')
    BONE_KW = ['ossification', 'complement', 'osteo', 'bone', 'inflamm', 'cytokine']
    for r in S.run("\n        MATCH (dem:Disease)-[:ASSOCIATES_DaG]-(g2:Gene)\n              -[:PARTICIPATES_GpPW]-(pw:Pathway)-[:PARTICIPATES_GpPW]-(spp1:Gene {identifier:$g})\n        WHERE (toLower(dem.name) CONTAINS 'dementia' OR toLower(dem.name) CONTAINS 'alzheimer')\n          AND g2.name IN $keep\n        RETURN DISTINCT dem.identifier AS dem_id, dem.name AS dem_name,\n                        g2.identifier AS g2_id, g2.name AS g2_name,\n                        pw.identifier AS pw_id, pw.name AS pw_name\n    ", g=SPP1_ID, keep=list(KEEP_DEM_PARTNERS)).data():
        pwl = r['pw_name'].lower()
        if not any((k in pwl for k in BONE_KW)):
            continue
        dem = add_node('Disease', r['dem_id'], r['dem_name'], role='Dementia')
        g2 = add_node('Gene', r['g2_id'], r['g2_name'], role='inflammatory')
        pw = add_node('Pathway', r['pw_id'], r['pw_name'])
        add_edge(dem, g2, 'ASSOCIATES_DaG', 'DaG')
        add_edge(g2, pw, 'PARTICIPATES_GpPW', 'GpPW')
        add_edge(pw, f'Gene:{SPP1_ID}', 'PARTICIPATES_GpPW', 'GpPW')
driver.close()
nodes_df = pd.DataFrame(list(nodes.values())).drop_duplicates(subset=['id']).reset_index(drop=True)
edges_df = pd.DataFrame(list(edges.values())).drop_duplicates().reset_index(drop=True)
used = set(edges_df['source']) | set(edges_df['target'])
nodes_df = nodes_df[nodes_df['id'].isin(used)].reset_index(drop=True)
NODES_PATH = os.path.join(OUT_DIR, 'spp1_pm_osteoporosis_nodes.csv')
EDGES_PATH = os.path.join(OUT_DIR, 'spp1_pm_osteoporosis_edges.csv')
nodes_df.to_csv(NODES_PATH, index=False)
edges_df.to_csv(EDGES_PATH, index=False)
print(f'Saved nodes: {NODES_PATH} (n={len(nodes_df)})')
print(f'Saved edges: {EDGES_PATH} (n={len(edges_df)})')
print('Node type counts:')
print(nodes_df['type'].value_counts().to_string())
