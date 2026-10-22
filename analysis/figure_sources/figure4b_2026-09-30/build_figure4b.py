# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Prepare Figure 4B from the saved 22 April 2026 Cytoscape export.

Keep every direct insomnia/paresthesia -> gene -> dementia/Alzheimer's
route in that export. Keep disease identifiers distinct. No live queries.
"""
import os
from pathlib import Path
import csv, json, collections, random, math, hashlib
OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent
nodes = list(csv.DictReader((OUT / 'source_cytoscape_nodes.csv').open()))
edges = list(csv.DictReader((OUT / 'source_cytoscape_edges.csv').open()))
lookup = {n['id']: n for n in nodes}
rows = list(csv.DictReader((ROOT / 'audit_2026-09-30/excel_stage/Q3_Symptom_Gene_Dementia.csv').open()))
I = 'Symptom:D007319'
P = 'Symptom:D010292'
sg = [e for e in edges if e['source'] in {I, P} and e['interaction'] == 'ASSOCIATES_GaS']
geneids = {e['target'] for e in sg}
diseaseids = {r['disease_id'] for r in rows}
gd = [e for e in edges if e['source'] in geneids and e['target'] in diseaseids and (e['interaction'] == 'ASSOCIATES_DaG') and (e['path_class'] == 'Gene-Disease')]
selected = sg + gd
assert len(geneids) == 15 and len(diseaseids) == 8 and (len(sg) == 17) and (len(gd) == 33)
expected = {(r['symptom_id'], r['gene_id'], r['disease_id']) for r in rows}
actual = {(a['source'], a['target'], b['target']) for a in sg for b in gd if a['target'] == b['source']}
assert expected == actual and len(actual) == 36
membership = {g: {e['source'] for e in sg if e['target'] == g} for g in geneids}
groups = [[g for g in sorted(geneids) if membership[g] == s] for s in [{I}, {I, P}, {P}]]
assert list(map(len, groups)) == [7, 2, 6]
reg = [(e['source'], e['target']) for e in gd]
names = {i: n['name'] for i, n in lookup.items()}

def score(gs, ds):
    a = {g: i for i, g in enumerate(sum(gs, []))}
    b = {d: i for i, d in enumerate(ds)}
    crosses = sum(((a[u] - a[x]) * (b[v] - b[y]) < 0 for i, (u, v) in enumerate(reg) for x, y in reg[i + 1:]))
    span = sum((abs(a[u] / 14 - b[v] / 7) for u, v in reg))
    return crosses * 5 + span
rng = random.Random(42)
targets = sorted(diseaseids)
best = (score(groups, targets), [g[:] for g in groups], targets[:])
for restart in range(5):
    gs = [g[:] for g in groups]
    ds = targets[:]
    if restart:
        for g in gs:
            rng.shuffle(g)
        rng.shuffle(ds)
    sc = score(gs, ds)
    for j in range(4500):
        choices = gs + [ds, ds]
        changed = rng.choice(choices)
        a, b = rng.sample(range(len(changed)), 2)
        changed[a], changed[b] = (changed[b], changed[a])
        new = score(gs, ds)
        temp = max(0.06, 10 * (1 - j / 4500) ** 2)
        if new < sc or rng.random() < math.exp(min(0, (sc - new) / temp)):
            sc = new
            if sc < best[0]:
                best = (sc, [g[:] for g in gs], ds[:])
        else:
            changed[a], changed[b] = (changed[b], changed[a])
order = sum(best[1], [])
target_order = best[2]
short = {'Disease:DOID:10652': "Alzheimer's\ndisease", 'Disease:DOID:1307': 'Dementia', 'Disease:DOID:12217': 'Lewy body\ndementia', 'Disease:DOID:9255': 'Frontotemporal\ndementia', 'Disease:DOID:8725': 'Vascular\ndementia', 'Disease:DOID:0060213': 'FTD/ALS 1', 'Disease:DOID:0111246': 'ALS-PDC 1', 'Disease:DOID:0050881': 'IBMPFD'}
net = {'data': {'name': 'Figure 4B | Complete direct symptom routes | 15 genes | April 2026'}, 'elements': {'nodes': [], 'edges': []}}
for gid in order:
    n = lookup[gid]
    net['elements']['nodes'].append({'data': {'id': gid, 'name': n['name'], 'label': n['name'], 'node_type': 'Gene', 'role': 'bridge_gene', 'source_node_id': gid, 'symptoms': '; '.join((names[s] for s in sorted(membership[gid])))}, 'position': {'x': 330, 'y': 45 + order.index(gid) * 40}})
for sid, y in [(I, 185), (P, 485)]:
    count = sum((sid in membership[g] for g in geneids))
    net['elements']['nodes'].append({'data': {'id': sid, 'name': lookup[sid]['name'], 'label': lookup[sid]['name'], 'node_type': 'Symptom', 'role': 'symptom_of_interest', 'source_node_id': sid}, 'position': {'x': 15, 'y': y}})
for i, did in enumerate(target_order):
    net['elements']['nodes'].append({'data': {'id': did, 'name': lookup[did]['name'], 'label': short[did], 'node_type': 'Disease', 'role': 'disease_of_interest', 'source_node_id': did}, 'position': {'x': 635, 'y': 45 + i * 80}})
for i, e in enumerate(selected):
    net['elements']['edges'].append({'data': {'id': f'e{i + 1:03d}', 'source': e['source'], 'target': e['target'], 'predicate': e['interaction'], 'interaction': e['interaction'], 'path_class': e['path_class'], 'source_record': 'Saved April 22 2026 CSV export'}})
defs = {'NETWORK_BACKGROUND_PAINT': '#FFFFFF', 'NODE_SHAPE': 'ELLIPSE', 'NODE_WIDTH': 30, 'NODE_HEIGHT': 30, 'NODE_FILL_COLOR': '#CEE4EF', 'NODE_BORDER_PAINT': '#4E7287', 'NODE_BORDER_WIDTH': 1.2, 'NODE_LABEL_COLOR': '#162E3E', 'NODE_LABEL_FONT_SIZE': 30, 'NODE_LABEL_FONT_FACE': 'Arial-ItalicMT,plain,30', 'NODE_LABEL_WIDTH': 230, 'NODE_LABEL_POSITION': 'W,E,r,-9.00,0.00', 'NODE_LABEL_BACKGROUND_COLOR': '#FFFFFF', 'NODE_LABEL_BACKGROUND_SHAPE': 'ROUND_RECTANGLE', 'NODE_LABEL_BACKGROUND_TRANSPARENCY': 255, 'NODE_VISIBLE': True, 'NODE_TRANSPARENCY': 255, 'EDGE_WIDTH': 1.6, 'EDGE_SOURCE_ARROW_SHAPE': 'NONE', 'EDGE_TARGET_ARROW_SHAPE': 'NONE', 'EDGE_TRANSPARENCY': 220, 'EDGE_VISIBLE': True, 'EDGE_LINE_TYPE': 'SOLID', 'EDGE_LABEL': '', 'EDGE_CURVED': False}
style = {'title': 'Figure4B_Compact_15genes_20260930', 'defaults': [{'visualProperty': k, 'value': v} for k, v in defs.items()], 'mappings': [{'mappingType': 'passthrough', 'mappingColumn': 'label', 'mappingColumnType': 'String', 'visualProperty': 'NODE_LABEL'}]}
for vp in ['EDGE_STROKE_UNSELECTED_PAINT', 'EDGE_UNSELECTED_PAINT']:
    style['mappings'].append({'mappingType': 'discrete', 'mappingColumn': 'predicate', 'mappingColumnType': 'String', 'visualProperty': vp, 'map': [{'key': 'ASSOCIATES_GaS', 'value': '#347DA5'}, {'key': 'ASSOCIATES_DaG', 'value': '#B2BBC1'}]})
annotations = []
for i, (text, color) in enumerate([('Symptom–gene association (GaS)', '#347DA5'), ('Gene–disease association (DaG)', '#A6B1B9')]):
    y = 672 + i * 44
    annotations.extend([{'kind': 'shape', 'newName': f'Edge key {i}', 'type': 'RECTANGLE', 'x': '40', 'y': str(y + 16), 'width': '55', 'height': '4', 'fillColor': color, 'borderThickness': '0', 'opacity': '100', 'canvas': 'foreground'}, {'kind': 'text', 'newName': f'Edge label {i}', 'text': text, 'x': '113', 'y': str(y), 'fontSize': '28', 'fontFamily': 'Arial', 'fontStyle': 'plain', 'color': '#263B48', 'canvas': 'foreground'}])
for f, data in [('network.json', net), ('style.json', style), ('annotations.json', annotations)]:
    (OUT / f).write_text(json.dumps(data, indent=2) + '\n')
with (OUT / 'Figure4B_nodes.csv').open('w', newline='') as f:
    fields = ['node_id', 'name', 'display_label', 'node_type']
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(({'node_id': n['data']['id'], 'name': n['data']['name'], 'display_label': n['data']['label'].replace('\n', ' / '), 'node_type': n['data']['node_type']} for n in net['elements']['nodes']))
with (OUT / 'Figure4B_edges.csv').open('w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(selected[0]))
    w.writeheader()
    w.writerows(selected)
with (OUT / 'Figure4B_complete_routes.csv').open('w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
caption = "b, Symptom–gene–disease association network for insomnia and paresthesia, showing all 15 intermediate genes and eight distinct dementia-related disease nodes recovered from the saved 22 April 2026 SPOKE export. Insomnia and paresthesia connect through 9 and 8 genes, respectively, with HLA-DRB1 and PRNP shared between them. Pink, light-blue and red circles indicate symptoms, genes and diseases. Blue lines denote symptom–gene associations (GaS), and gray lines gene–disease associations (DaG). Gene symbols are italicized. No gene whitelist was applied; indirect routes through musculoskeletal pain are outside this panel's scope. FTD/ALS 1, frontotemporal dementia and/or amyotrophic lateral sclerosis 1; ALS-PDC 1, amyotrophic lateral sclerosis–parkinsonism/dementia complex 1; IBMPFD, inclusion body myopathy with Paget disease of bone and frontotemporal dementia."
(OUT / 'Figure4B_caption.txt').write_text(caption + '\n')
provenance = {'data_version': 'Saved April 22 2026 SPOKE export; no new queries', 'nodes': 25, 'genes': 15, 'symptoms': 2, 'distinct_disease_endpoints': 8, 'edges': 50, 'symptom_gene_edges': 17, 'unique_gene_disease_edges': 33, 'complete_routes': 36, 'insomnia_genes': 9, 'insomnia_associations': 16, 'paresthesia_genes': 8, 'paresthesia_associations': 20, 'shared_genes': ['HLA-DRB1', 'PRNP'], 'scope': 'All direct insomnia/paresthesia–gene–dementia/Alzheimer disease routes; indirect musculoskeletal-pain routes excluded', 'comments_addressed': ['0: italic gene symbols', '235: remove manual gene whitelist and explain distinct-gene vs association counts'], 'layout': 'Three layers, gene ordering constrained by symptom membership; seeded crossing reduction. Coordinates carry no biological-distance meaning.', 'source_sha256': {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in OUT.glob('source_*.csv')}}
(OUT / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
print(json.dumps({'nodes': 25, 'edges': 50, 'genes': 15, 'gene_order': [names[g] for g in order], 'disease_order': [names[d] for d in target_order]}, indent=2))
