# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os

def extract_intermediates_from_paths(G, source_node, target_node, metapath, max_paths=200):
    """
    Walk paths following exact metapath edge sequence (MultiGraph-safe).
    Extracts ALL intermediate nodes, not just genes.
    """

    def walk(current_node, step, path):
        if step == len(metapath):
            if current_node == target_node:
                results.append(list(path))
            return
        if len(results) >= max_paths:
            return
        edge_type = metapath[step]
        for neighbor in G.neighbors(current_node):
            if len(results) >= max_paths:
                return
            edge_dict = G[current_node][neighbor]
            has_match = any((edata.get('predicate') == edge_type for edata in edge_dict.values()))
            if has_match:
                path.append(neighbor)
                walk(neighbor, step + 1, path)
                path.pop()
    results = []
    walk(source_node, 0, [source_node])
    intermediates = []
    for path in results:
        for node in path[1:-1]:
            intermediates.append({'node_id': node, 'node_name': G.nodes[node].get('name', node), 'node_type': G.nodes[node].get('node_type', 'Unknown'), 'metapath': ' -> '.join(metapath)})
    return pd.DataFrame(intermediates)
TOP_K = 5
top_metapaths = dwpc_nonzero.sort_values('dwpc', ascending=False).groupby('target_name').head(TOP_K)
print(f'Extracting intermediates from top {TOP_K} metapaths per target...')
print(f'Total metapaths to process: {len(top_metapaths)}')
all_intermediates = []
for _, row in tqdm(top_metapaths.iterrows(), total=len(top_metapaths)):
    metapath = row['metapath'].split(' → ')
    inter_df = extract_intermediates_from_paths(G, source_node, row['target_id'], metapath, max_paths=200)
    if len(inter_df) > 0:
        inter_df['target_name'] = row['target_name']
        inter_df['dwpc'] = row['dwpc']
        all_intermediates.append(inter_df)
intermediates_df = pd.concat(all_intermediates, ignore_index=True) if all_intermediates else pd.DataFrame()
print(f'\nTotal intermediate node occurrences: {len(intermediates_df)}')
if len(intermediates_df) == 0:
    print('No intermediates found — check metapath walking logic.')
else:
    print('\nIntermediate node types found:')
    print(intermediates_df['node_type'].value_counts().to_string())
    genes_df = intermediates_df[intermediates_df['node_type'] == 'Gene'].copy()
    print(f'\nGene intermediates: {len(genes_df)}')
    non_genes = intermediates_df[intermediates_df['node_type'] != 'Gene']
    print(f'Non-gene intermediates: {len(non_genes)}')
    for nt in non_genes['node_type'].unique():
        subset = non_genes[non_genes['node_type'] == nt]
        print(f"\n  {nt} ({len(subset['node_name'].unique())} unique):")
        print(f"    {', '.join(subset['node_name'].unique()[:10])}")
    if len(genes_df) > 0:
        gene_summary = genes_df.groupby(['node_id', 'node_name']).agg(n_targets=('target_name', 'nunique'), n_metapaths=('metapath', 'nunique'), total_occurrences=('node_name', 'count'), targets=('target_name', lambda x: ', '.join(sorted(x.unique()))), max_dwpc=('dwpc', 'max')).reset_index().rename(columns={'node_name': 'gene_name'}).sort_values(['n_targets', 'total_occurrences'], ascending=False)
        print(f'\nUnique genes found: {len(gene_summary)}')
        print(f"Genes appearing in ALL 3 target paths: {len(gene_summary[gene_summary['n_targets'] == 3])}")
        print(f"Genes appearing in 2+ target paths: {len(gene_summary[gene_summary['n_targets'] >= 2])}")
        print('\nTop 30 recurring genes:')
        print(gene_summary.head(30).to_string(index=False))
        gene_summary.to_csv('dwpc_recurring_genes.csv', index=False)
    else:
        print('\nNo Gene intermediates in top metapaths.')
        print('Top metapaths route through Disease/Compound/Anatomy nodes instead.')
        print('Consider increasing TOP_K or filtering for gene-containing metapaths.')
