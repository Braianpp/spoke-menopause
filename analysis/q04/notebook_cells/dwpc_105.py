# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
TOP_N = 200
pm_top = set(pm_genes.head(TOP_N)['gene_name'])
osteo_top = set(osteo_genes.head(TOP_N)['gene_name'])
dementia_top = set(dementia_genes.head(TOP_N)['gene_name'])
all_three = pm_top & osteo_top & dementia_top
pm_osteo = (pm_top & osteo_top) - dementia_top
pm_dementia = (pm_top & dementia_top) - osteo_top
osteo_dementia = (osteo_top & dementia_top) - pm_top
bone_brain = osteo_top & dementia_top
print(f'Top {TOP_N} genes per condition:\n')
print(f'All three (PM + Osteoporosis + Dementia): {len(all_three)}')
print(f'PM + Osteoporosis only: {len(pm_osteo)}')
print(f'PM + Dementia only: {len(pm_dementia)}')
print(f'Osteoporosis + Dementia only: {len(osteo_dementia)}')
print(f'Bone-Brain overlap (Osteo + Dementia, incl. triple): {len(bone_brain)}')
overlap_genes = bone_brain | pm_osteo | pm_dementia | all_three
overlap_data = []
for gene in overlap_genes:
    pm_rank = pm_genes.loc[pm_genes['gene_name'] == gene, 'rank'].values
    osteo_rank = osteo_genes.loc[osteo_genes['gene_name'] == gene, 'rank'].values
    dem_rank = dementia_genes.loc[dementia_genes['gene_name'] == gene, 'rank'].values
    pm_r = pm_rank[0] if len(pm_rank) > 0 else None
    osteo_r = osteo_rank[0] if len(osteo_rank) > 0 else None
    dem_r = dem_rank[0] if len(dem_rank) > 0 else None
    in_pm = gene in pm_top
    in_osteo = gene in osteo_top
    in_dementia = gene in dementia_top
    category = 'Other'
    if in_pm and in_osteo and in_dementia:
        category = 'All three'
    elif in_osteo and in_dementia:
        category = 'Bone-Brain (Osteo+Dementia)'
    elif in_pm and in_osteo:
        category = 'PM+Osteoporosis'
    elif in_pm and in_dementia:
        category = 'PM+Dementia'
    overlap_data.append({'gene_name': gene, 'pm_rank': pm_r, 'osteo_rank': osteo_r, 'dementia_rank': dem_r, 'category': category})
overlap_df = pd.DataFrame(overlap_data)
overlap_df['avg_rank'] = overlap_df[['pm_rank', 'osteo_rank', 'dementia_rank']].mean(axis=1)
cat_order = {'All three': 0, 'Bone-Brain (Osteo+Dementia)': 1, 'PM+Osteoporosis': 2, 'PM+Dementia': 3}
overlap_df['cat_sort'] = overlap_df['category'].map(cat_order)
overlap_df = overlap_df.sort_values(['cat_sort', 'avg_rank']).drop(columns='cat_sort')
print(f"\n{'=' * 70}")
print(f'TRIPLE OVERLAP GENES (PM + Osteoporosis + Dementia): {len(all_three)}')
print(f"{'=' * 70}")
print(overlap_df[overlap_df['category'] == 'All three'][['gene_name', 'pm_rank', 'osteo_rank', 'dementia_rank']].to_string(index=False))
print(f"\n{'=' * 70}")
print(f'BONE-BRAIN GENES (Osteoporosis + Dementia): {len(osteo_dementia)}')
print(f"{'=' * 70}")
print(overlap_df[overlap_df['category'] == 'Bone-Brain (Osteo+Dementia)'].head(20)[['gene_name', 'pm_rank', 'osteo_rank', 'dementia_rank']].to_string(index=False))
