# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""
Generate Cytoscape-ready node and edge tables for a focused
Symptom–Gene–Disease subnetwork linking selected perimenopausal symptoms
to dementia via SPOKE (Neo4j).

Outputs (CSV):
- cytoscape_nodes.csv
- cytoscape_edges.csv

Run on leo:
  python3 ${SPOKE_DATA_ROOT}/biorouter/scripts/q_cytoscape_symptom_gene_dementia_v1.py
"""
import os
import os
import sys
import pandas as pd
from neo4j import GraphDatabase
SPOKE_DEV_URI = os.environ.get('SPOKE_DEV_URI')
SPOKE_User = os.environ.get('SPOKE_User')
SPOKE_Password = os.environ.get('SPOKE_Password')
_candidates = [os.path.dirname(__file__), os.path.dirname(os.path.dirname(__file__)), '${SPOKE_DATA_ROOT}', '${SPOKE_DATA_ROOT}/biorouter']
for p in _candidates:
    if p and os.path.isdir(p) and (p not in sys.path):
        sys.path.append(p)
try:
    from spoke_utils import SPOKE_DEV_URI as _URI, SPOKE_User as _USER, SPOKE_Password as _PWD
    SPOKE_DEV_URI = SPOKE_DEV_URI or _URI
    SPOKE_User = SPOKE_User or _USER
    SPOKE_Password = SPOKE_Password or _PWD
except Exception:
    pass
if not (SPOKE_DEV_URI and SPOKE_User and SPOKE_Password):
    print('Missing SPOKE connection settings; set SPOKE_DEV_URI, SPOKE_User, SPOKE_Password.', file=sys.stderr)
    sys.exit(1)
SYMPTOMS = {'Insomnia': 'D007319', 'Paresthesia': 'D010292', 'Musculoskeletal Pain': 'D059352'}
symptom_ids = [SYMPTOMS['Insomnia'], SYMPTOMS['Paresthesia']]
msk_id = SYMPTOMS['Musculoskeletal Pain']
OUT_DIR = os.environ.get('CYTO_OUT', '${SPOKE_DATA_ROOT}/biorouter/exports')
os.makedirs(OUT_DIR, exist_ok=True)
driver = GraphDatabase.driver(SPOKE_DEV_URI, auth=(SPOKE_User, SPOKE_Password))
CYPHER_DIRECT = {'query': "\n        MATCH (s:Symptom)-[:ASSOCIATES_GaS]-(g:Gene)-[:ASSOCIATES_DaG]-(d:Disease)\n        WHERE s.identifier IN $symptom_ids\n          AND (toLower(d.name) CONTAINS 'dementia' OR toLower(d.name) CONTAINS 'alzheimer')\n        RETURN DISTINCT\n          s.identifier AS symptom_id,\n          s.name       AS symptom_name,\n          g.identifier AS gene_id,\n          g.name       AS gene_name,\n          d.identifier AS disease_id,\n          d.name       AS disease_name\n        ", 'params': {'symptom_ids': symptom_ids}}
CYPHER_INDIRECT = {'query': "\n        MATCH (d2:Disease)-[:PRESENTS_DpS]->(s2:Symptom)\n        WHERE s2.identifier = $msk_id\n        WITH d2, s2\n        MATCH (d2)-[:ASSOCIATES_DaG]-(g:Gene)-[:ASSOCIATES_DaG]-(d1:Disease)\n        WHERE (toLower(d1.name) CONTAINS 'dementia' OR toLower(d1.name) CONTAINS 'alzheimer')\n        RETURN DISTINCT\n          s2.identifier AS symptom_id,\n          s2.name       AS symptom_name,\n          d2.identifier AS presenting_disease_id,\n          d2.name       AS presenting_disease_name,\n          g.identifier  AS gene_id,\n          g.name        AS gene_name,\n          d1.identifier AS disease_id,\n          d1.name       AS disease_name\n        ", 'params': {'msk_id': msk_id}}
with driver.session() as session:
    direct_df = pd.DataFrame(session.run(CYPHER_DIRECT['query'], **CYPHER_DIRECT['params']).data())
    indirect_df = pd.DataFrame(session.run(CYPHER_INDIRECT['query'], **CYPHER_INDIRECT['params']).data())
direct_df = direct_df if not direct_df.empty else pd.DataFrame(columns=['symptom_id', 'symptom_name', 'gene_id', 'gene_name', 'disease_id', 'disease_name'])
indirect_df = indirect_df if not indirect_df.empty else pd.DataFrame(columns=['symptom_id', 'symptom_name', 'presenting_disease_id', 'presenting_disease_name', 'gene_id', 'gene_name', 'disease_id', 'disease_name'])
nodes = []
for _, row in pd.concat([direct_df[['symptom_id', 'symptom_name']], indirect_df[['symptom_id', 'symptom_name']]], ignore_index=True).drop_duplicates().iterrows():
    nodes.append({'id': f'Symptom:{row.symptom_id}', 'name': row.symptom_name, 'type': 'Symptom', 'role': 'symptom_of_interest'})
if not indirect_df.empty:
    for _, row in indirect_df[['presenting_disease_id', 'presenting_disease_name']].drop_duplicates().iterrows():
        nodes.append({'id': f'Disease:{row.presenting_disease_id}', 'name': row.presenting_disease_name, 'type': 'Disease', 'role': 'intermediate_disease'})
for _, row in pd.concat([direct_df[['disease_id', 'disease_name']], indirect_df[['disease_id', 'disease_name']]], ignore_index=True).drop_duplicates().iterrows():
    nodes.append({'id': f'Disease:{row.disease_id}', 'name': row.disease_name, 'type': 'Disease', 'role': 'disease_of_interest'})
for _, row in pd.concat([direct_df[['gene_id', 'gene_name']], indirect_df[['gene_id', 'gene_name']]], ignore_index=True).drop_duplicates().iterrows():
    nodes.append({'id': f'Gene:{row.gene_id}', 'name': row.gene_name, 'type': 'Gene', 'role': 'bridge_gene'})
nodes_df = pd.DataFrame(nodes).drop_duplicates(subset=['id']).reset_index(drop=True)
edges = []
if not direct_df.empty:
    for _, row in direct_df.iterrows():
        edges.append({'source': f'Symptom:{row.symptom_id}', 'target': f'Gene:{row.gene_id}', 'interaction': 'ASSOCIATES_GaS', 'path_class': 'Symptom-Gene'})
if not indirect_df.empty:
    for _, row in indirect_df.iterrows():
        edges.append({'source': f'Disease:{row.presenting_disease_id}', 'target': f'Symptom:{row.symptom_id}', 'interaction': 'PRESENTS_DpS', 'path_class': 'Disease-Symptom'})
        edges.append({'source': f'Disease:{row.presenting_disease_id}', 'target': f'Gene:{row.gene_id}', 'interaction': 'ASSOCIATES_DaG', 'path_class': 'Disease-Gene'})
for _, row in pd.concat([direct_df[['gene_id', 'disease_id']], indirect_df[['gene_id', 'disease_id']]], ignore_index=True).drop_duplicates().iterrows():
    edges.append({'source': f'Gene:{row.gene_id}', 'target': f'Disease:{row.disease_id}', 'interaction': 'ASSOCIATES_DaG', 'path_class': 'Gene-Disease'})
edges_df = pd.DataFrame(edges).drop_duplicates().reset_index(drop=True)
nodes_path = os.path.join(OUT_DIR, 'cytoscape_nodes.csv')
edges_path = os.path.join(OUT_DIR, 'cytoscape_edges.csv')
nodes_df.to_csv(nodes_path, index=False)
edges_df.to_csv(edges_path, index=False)
print(f'Saved nodes to: {nodes_path}  (n={len(nodes_df)})')
print(f'Saved edges to: {edges_path}  (n={len(edges_df)})')
