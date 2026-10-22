# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib_venn import venn3
pm_emb = pd.read_csv('${SPOKE_DATA_ROOT}/premature_menopause_embedding.csv')
osteo_emb = pd.read_csv('${SPOKE_DATA_ROOT}/osteoporosis_embedding.csv')
dementia_emb = pd.read_csv('${SPOKE_DATA_ROOT}/dementia_embedding.csv')
print(f'PM embedding: {len(pm_emb)} nodes')
print(f'Osteoporosis embedding: {len(osteo_emb)} nodes')
print(f'Dementia embedding: {len(dementia_emb)} nodes')

def get_gene_rankings(emb_df, G):
    """Filter embedding to gene nodes and rank by PageRank score."""
    gene_rows = []
    for _, row in emb_df.iterrows():
        node_raw = str(row['node'])
        if node_raw.isdigit():
            node = f'Gene:{node_raw}'
        else:
            node = node_raw
        if node in G.nodes and G.nodes[node].get('node_type') == 'Gene':
            gene_rows.append({'node': node, 'gene_name': G.nodes[node].get('name', node), 'pagerank': row['pagerank']})
    gene_df = pd.DataFrame(gene_rows).sort_values('pagerank', ascending=False).reset_index(drop=True)
    gene_df['rank'] = range(1, len(gene_df) + 1)
    return gene_df
pm_genes = get_gene_rankings(pm_emb, G)
osteo_genes = get_gene_rankings(osteo_emb, G)
dementia_genes = get_gene_rankings(dementia_emb, G)
print(f'\nGenes ranked - PM: {len(pm_genes)}, Osteoporosis: {len(osteo_genes)}, Dementia: {len(dementia_genes)}')
print(f'\nTop 10 genes by PageRank:')
print(f'\n--- Premature Menopause ---')
print(pm_genes.head(10)[['rank', 'gene_name', 'pagerank']].to_string(index=False))
print(f'\n--- Osteoporosis ---')
print(osteo_genes.head(10)[['rank', 'gene_name', 'pagerank']].to_string(index=False))
print(f'\n--- Dementia ---')
print(dementia_genes.head(10)[['rank', 'gene_name', 'pagerank']].to_string(index=False))
