# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
calculator = NetworkDWPC(G)
source_node = 'Disease:DOID:10787'
source_type = G.nodes[source_node]['node_type']
source_name = G.nodes[source_node].get('name', source_node)
targets = {'Disease:DOID:8725': 'Vascular dementia', 'Disease:DOID:9255': 'Frontotemporal dementia', 'Disease:DOID:14330': "Parkinson's disease"}
all_results = []
for target_node, target_label in tqdm(targets.items()):
    if target_node not in G.nodes:
        print(f'WARNING: {target_label} ({target_node}) not found in G_small, skipping.')
        continue
    target_type = G.nodes[target_node]['node_type']
    metapaths = calculator.find_metapaths(source_type, target_type, max_length=3)
    for metapath in metapaths:
        dwpc = calculator.compute_dwpc(source_node, target_node, metapath)
        path_count = calculator.compute_simple_path_count(source_node, target_node, metapath)
        all_results.append({'source_id': source_node, 'source_name': 'Premature menopause', 'target_id': target_node, 'target_name': target_label, 'metapath': ' → '.join(metapath), 'length': len(metapath), 'dwpc': dwpc, 'path_count': path_count})
dwpc_df = pd.DataFrame(all_results)
print(f'Total metapath results: {len(dwpc_df)}')
dwpc_df.sort_values('dwpc', ascending=False).head(20)
