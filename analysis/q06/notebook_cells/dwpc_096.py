# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
calculator = NetworkDWPC(G)
pairs = [('Disease:DOID:10787', 'Premature menopause', 'Disease:DOID:8645', 'Subacute delirium'), ('Disease:DOID:10787', 'Premature menopause', 'Disease:DOID:1307', 'Dementia'), ('Disease:DOID:8645', 'Subacute delirium', 'Disease:DOID:1307', 'Dementia')]
all_results = []
for source_id, source_label, target_id, target_label in pairs:
    if source_id not in G.nodes:
        print(f'WARNING: {source_label} ({source_id}) not found, skipping.')
        continue
    if target_id not in G.nodes:
        print(f'WARNING: {target_label} ({target_id}) not found, skipping.')
        continue
    source_type = G.nodes[source_id]['node_type']
    target_type = G.nodes[target_id]['node_type']
    metapaths = calculator.find_metapaths(source_type, target_type, max_length=3)
    pair_label = f'{source_label} → {target_label}'
    print(f'Computing DWPC for {pair_label} ({len(metapaths)} metapaths)...')
    for metapath in tqdm(metapaths, desc=pair_label):
        dwpc = calculator.compute_dwpc(source_id, target_id, metapath)
        path_count = calculator.compute_simple_path_count(source_id, target_id, metapath)
        all_results.append({'source_id': source_id, 'source_name': source_label, 'target_id': target_id, 'target_name': target_label, 'pair': pair_label, 'metapath': ' → '.join(metapath), 'length': len(metapath), 'dwpc': dwpc, 'path_count': path_count})
dwpc_del_df = pd.DataFrame(all_results)
dwpc_del_nonzero = dwpc_del_df[dwpc_del_df['dwpc'] > 0].copy()
print(f'\nTotal metapath results: {len(dwpc_del_df)}')
print(f'Non-zero: {len(dwpc_del_nonzero)}')
for pair in dwpc_del_df['pair'].unique():
    subset = dwpc_del_nonzero[dwpc_del_nonzero['pair'] == pair]
    print(f"\n{'=' * 60}")
    print(f"{pair}: {len(subset)} active metapaths, cumulative DWPC = {subset['dwpc'].sum():.4f}")
    print(f"{'=' * 60}")
    print(subset.sort_values('dwpc', ascending=False).head(10)[['metapath', 'dwpc', 'path_count']].to_string(index=False))
