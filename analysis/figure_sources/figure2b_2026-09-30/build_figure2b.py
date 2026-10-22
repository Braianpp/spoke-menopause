# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Prepare Figure 2B for native Cytoscape rendering from audited source tables.

No graph loading or new biological queries. Coordinates reduce edge crossings;
they do not represent biological distances or computed functional groups.
"""
import os
from pathlib import Path
import csv
import json
import random
import math
import hashlib
OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent
STAGE = ROOT / 'audit_2026-09-30/excel_stage'
nodes = list(csv.DictReader((STAGE / 'Q1_Fig2b_Network_Nodes.csv').open()))
edges = list(csv.DictReader((STAGE / 'Q1_Fig2b_Network_Edges.csv').open()))
PM = 'Disease:DOID:10787'
DEM = 'Disease:DOID:1307'
names = {n['node_id']: n['name'] for n in nodes}
left = sorted({e['target'] for e in edges if e['source'] == PM}, key=names.get)
right = sorted({e['source'] for e in edges if e['target'] == DEM}, key=names.get)
reg = [(e['source'], e['target']) for e in edges if e['predicate'] != 'ASSOCIATES_DaG']
assert len(nodes) == 36 and len(edges) == 67
assert len(left) + len(right) == 34 and (not set(left) & set(right))
assert set(left) | set(right) | {PM, DEM} == set(names)

def cost(a, b):
    ai = {n: i for i, n in enumerate(a)}
    bi = {n: i for i, n in enumerate(b)}
    crossing = sum(((ai[u] - ai[x]) * (bi[v] - bi[y]) < 0 for i, (u, v) in enumerate(reg) for x, y in reg[i + 1:]))
    span = sum((abs(ai[u] / (len(a) - 1) - bi[v] / (len(b) - 1)) for u, v in reg))
    return crossing * 10 + span
rng = random.Random(42)
best = (cost(left, right), left[:], right[:])
for restart in range(12):
    a = left[:]
    b = right[:]
    if restart:
        rng.shuffle(a)
        rng.shuffle(b)
    score = cost(a, b)
    for i in range(6000):
        target = a if rng.random() < 0.5 else b
        p, q = rng.sample(range(len(target)), 2)
        target[p], target[q] = (target[q], target[p])
        proposal = cost(a, b)
        temperature = max(0.08, 16 * (1 - i / 6000) ** 2)
        if proposal < score or rng.random() < math.exp(min(0, (score - proposal) / temperature)):
            score = proposal
            if score < best[0]:
                best = (score, a[:], b[:])
        else:
            target[p], target[q] = (target[q], target[p])
left, right = best[1:]
positions = {PM: (50, 430), DEM: (1080, 430)}
for group, x in [(left, 340), (right, 790)]:
    for i, n in enumerate(group):
        positions[n] = (x, 70 + i * 720 / (len(group) - 1))
network = {'data': {'name': 'Figure 2B | Complete regulation network | 34 genes'}, 'elements': {'nodes': [], 'edges': []}}
for n in nodes:
    nid = n['node_id']
    d = {'id': nid, 'name': n['name'], 'label': 'Premature\nmenopause' if nid == PM else 'Dementia' if nid == DEM else n['name'], 'node_type': n['node_type'], 'n_paths': int(n['n_paths']), 'source_node_id': nid, 'role': 'Premature menopause' if nid == PM else 'Dementia' if nid == DEM else 'First intermediate' if nid in left else 'Second intermediate'}
    x, y = positions[nid]
    network['elements']['nodes'].append({'data': d, 'position': {'x': x, 'y': y}})
for i, e in enumerate(edges):
    d = {'id': f'e{i + 1:03d}', 'source': e['source'], 'target': e['target'], 'interaction': e['predicate'], 'predicate': e['predicate'], 'n_paths': int(e['n_paths']), 'source_name': e['source_name'], 'target_name': e['target_name']}
    network['elements']['edges'].append({'data': d})
defaults = {'NETWORK_BACKGROUND_PAINT': '#FFFFFF', 'NODE_SHAPE': 'ROUND_RECTANGLE', 'NODE_WIDTH': 122.0, 'NODE_HEIGHT': 34.0, 'NODE_FILL_COLOR': '#CEE4EF', 'NODE_BORDER_PAINT': '#4E7287', 'NODE_BORDER_WIDTH': 1.2, 'NODE_LABEL_COLOR': '#162E3E', 'NODE_LABEL_FONT_SIZE': 19, 'NODE_LABEL_FONT_FACE': 'Arial-ItalicMT,plain,19', 'NODE_LABEL_WIDTH': 145.0, 'NODE_TRANSPARENCY': 255, 'NODE_VISIBLE': True, 'EDGE_WIDTH': 1.65, 'EDGE_STROKE_UNSELECTED_PAINT': '#B6BFC5', 'EDGE_UNSELECTED_PAINT': '#B6BFC5', 'EDGE_SOURCE_ARROW_SHAPE': 'NONE', 'EDGE_TARGET_ARROW_SHAPE': 'NONE', 'EDGE_LABEL': '', 'EDGE_TRANSPARENCY': 215, 'EDGE_VISIBLE': True, 'EDGE_CURVED': False}
style = {'title': 'Figure2B_34genes_20260930', 'defaults': [{'visualProperty': k, 'value': v} for k, v in defaults.items()], 'mappings': [{'mappingType': 'passthrough', 'mappingColumn': 'label', 'mappingColumnType': 'String', 'visualProperty': 'NODE_LABEL'}]}
colors = {'ASSOCIATES_DaG': '#BBC3C9', 'UPREGULATES_GPuG': '#C96532', 'UPREGULATES_KGuG': '#C96532', 'DOWNREGULATES_GPdG': '#347DA5', 'DOWNREGULATES_KGdG': '#347DA5'}
for vp, values in [('EDGE_STROKE_UNSELECTED_PAINT', colors), ('EDGE_UNSELECTED_PAINT', colors), ('EDGE_LINE_TYPE', {k: 'LONG_DASH' if '_KG' in k else 'SOLID' for k in colors}), ('EDGE_WIDTH', {k: 1.1 if k == 'ASSOCIATES_DaG' else 2.1 for k in colors})]:
    style['mappings'].append({'mappingType': 'discrete', 'mappingColumn': 'predicate', 'mappingColumnType': 'String', 'visualProperty': vp, 'map': [{'key': k, 'value': str(v)} for k, v in values.items()]})
for name, data in [('network.json', network), ('style.json', style), ('positions.json', positions)]:
    (OUT / name).write_text(json.dumps(data, indent=2))
for name in ['Q1_Fig2b_Network_Nodes.csv', 'Q1_Fig2b_Network_Edges.csv']:
    (OUT / name).write_bytes((STAGE / name).read_bytes())
manifest = {'source': 'Audited CSV exports from q1_pm_dementia_regulation_metapaths, 2026-06-11', 'nodes': 36, 'genes': 34, 'edges': 67, 'regulation_edges': len(reg), 'first_intermediates': len(left), 'second_intermediates': len(right), 'extraction': '4 regulation metapaths, maximum 50 paths per metapath; saved subset only', 'direction': 'Source/target are traversal order. Lines have no arrowheads because biological direction was not preserved separately.', 'layout': 'Two gene columns with deterministic crossing reduction, seed 42; coordinates have no biological-distance meaning.', 'source_sha256': {name: hashlib.sha256((STAGE / name).read_bytes()).hexdigest() for name in ['Q1_Fig2b_Network_Nodes.csv', 'Q1_Fig2b_Network_Edges.csv']}}
(OUT / 'provenance.json').write_text(json.dumps(manifest, indent=2))
print(json.dumps(manifest, indent=2))
