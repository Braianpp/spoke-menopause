# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os

def extract_genes_from_paths(G, source_node, target_node, metapath, max_paths=200):

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
            if any((edata.get('predicate') == edge_type for edata in edge_dict.values())):
                path.append(neighbor)
                walk(neighbor, step + 1, path)
                path.pop()
    results = []
    walk(source_node, 0, [source_node])
    genes = []
    for path in results:
        for node in path[1:-1]:
            if G.nodes[node].get('node_type') == 'Gene':
                genes.append({'node_id': node, 'gene_name': G.nodes[node].get('name', node), 'metapath': ' → '.join(metapath)})
    return pd.DataFrame(genes)
TOP_K = 5
top_metapaths = dwpc_del_nonzero.sort_values('dwpc', ascending=False).groupby('pair').head(TOP_K)
print(f'Extracting genes from top {TOP_K} metapaths per pair...')
all_genes = []
for _, row in tqdm(top_metapaths.iterrows(), total=len(top_metapaths)):
    metapath = row['metapath'].split(' → ')
    gene_df = extract_genes_from_paths(G, row['source_id'], row['target_id'], metapath, max_paths=200)
    if len(gene_df) > 0:
        gene_df['pair'] = row['pair']
        gene_df['source_name'] = row['source_name']
        gene_df['target_name'] = row['target_name']
        gene_df['dwpc'] = row['dwpc']
        all_genes.append(gene_df)
genes_del_df = pd.concat(all_genes, ignore_index=True) if all_genes else pd.DataFrame()
print(f'Total gene occurrences: {len(genes_del_df)}')
if len(genes_del_df) > 0:
    gene_pair_summary = genes_del_df.groupby(['node_id', 'gene_name']).agg(pairs=('pair', lambda x: set(x)), n_pairs=('pair', 'nunique'), total_occurrences=('gene_name', 'count')).reset_index().sort_values(['n_pairs', 'total_occurrences'], ascending=False)

    def label_overlap(pairs_set):
        in_pm_del = any(('Premature menopause → Subacute delirium' in p for p in pairs_set))
        in_pm_dem = any(('Premature menopause → Dementia' in p for p in pairs_set))
        in_del_dem = any(('Subacute delirium → Dementia' in p for p in pairs_set))
        if in_pm_del and in_pm_dem and in_del_dem:
            return 'All three pairs'
        elif in_pm_del and in_del_dem:
            return 'Menopause-Delirium + Delirium-Dementia'
        elif in_pm_dem and in_del_dem:
            return 'Menopause-Dementia + Delirium-Dementia'
        elif in_pm_del and in_pm_dem:
            return 'Menopause-Delirium + Menopause-Dementia'
        elif in_pm_del:
            return 'Menopause-Delirium only'
        elif in_pm_dem:
            return 'Menopause-Dementia only'
        else:
            return 'Delirium-Dementia only'
    gene_pair_summary['overlap'] = gene_pair_summary['pairs'].apply(label_overlap)
    print(f'\nGene overlap categories:')
    print(gene_pair_summary['overlap'].value_counts().to_string())
    print(f'\nGenes in ALL THREE pairs:')
    all3 = gene_pair_summary[gene_pair_summary['overlap'] == 'All three pairs']
    print(all3[['gene_name', 'total_occurrences']].to_string(index=False))
    print(f'\nGenes bridging Delirium (Meno-Del + Del-Dem):')
    bridge = gene_pair_summary[gene_pair_summary['overlap'] == 'Menopause-Delirium + Delirium-Dementia']
    print(bridge[['gene_name', 'total_occurrences']].head(20).to_string(index=False))
    print(f'\nTop 30 genes overall:')
    print(gene_pair_summary.head(30)[['gene_name', 'n_pairs', 'overlap', 'total_occurrences']].to_string(index=False))
    gene_pair_summary.to_csv('${SPOKE_DATA_ROOT}/dwpc_delirium_genes.csv', index=False)
