# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
import gseapy as gp
from matplotlib_venn import venn3
pm_del_genes = set(genes_del_df[genes_del_df['pair'].str.contains('Premature menopause → Subacute delirium')]['gene_name'])
pm_dem_genes = set(genes_del_df[genes_del_df['pair'].str.contains('Premature menopause → Dementia')]['gene_name'])
del_dem_genes = set(genes_del_df[genes_del_df['pair'].str.contains('Subacute delirium → Dementia')]['gene_name'])
triple = pm_del_genes & pm_dem_genes & del_dem_genes
delirium_bridge = pm_del_genes & del_dem_genes
print(f'PM→Delirium genes: {len(pm_del_genes)}')
print(f'PM→Dementia genes: {len(pm_dem_genes)}')
print(f'Delirium→Dementia genes: {len(del_dem_genes)}')
print(f'Triple overlap: {len(triple)}')
print(f'Delirium bridge genes (PM-Del + Del-Dem): {len(delirium_bridge)}')
enrich_genes = list(delirium_bridge) if len(delirium_bridge) >= 10 else list(pm_del_genes | del_dem_genes)
print(f'\nRunning enrichment on {len(enrich_genes)} genes')
gene_set_libraries = ['KEGG_2021_Human', 'Reactome_2022', 'GO_Biological_Process_2023', 'WikiPathway_2023_Human', 'Jensen_DISEASES']
enrichment_results = {}
for lib in gene_set_libraries:
    try:
        enr = gp.enrichr(gene_list=enrich_genes, gene_sets=lib, organism='human', outdir=None, cutoff=0.05)
        enrichment_results[lib] = enr.results[enr.results['Adjusted P-value'] < 0.05]
        print(f'  {lib}: {len(enrichment_results[lib])} significant terms')
    except Exception as e:
        print(f'  {lib}: FAILED - {e}')
all_enrichment_del = pd.concat([df.assign(database=db) for db, df in enrichment_results.items()], ignore_index=True)
axes_keywords = {'Cholinergic/Neurotransmitter': ['acetylchol', 'cholinerg', 'dopamin', 'serotonin', 'GABA', 'glutamat', 'neurotransmit', 'synapt'], 'Neuroinflammatory': ['inflam', 'cytokine', 'interleukin', 'TNF', 'NF-kB', 'chemokine', 'immune', 'toll-like', 'microglia', 'astrocyte'], 'Neurodegenerative': ['neuro', 'amyloid', 'tau', 'dementia', 'alzheimer', 'cognitive', 'brain', 'axon', 'myelin'], 'Hormonal/Endocrine': ['estrogen', 'steroid', 'hormone', 'gonad', 'ovari', 'androgen', 'aromatase', 'menopause'], 'Vascular/BBB': ['vascul', 'angiogen', 'endotheli', 'blood-brain', 'hypoxia', 'ischem', 'stroke', 'cerebrovascul']}

def classify_axis(term_name):
    term_lower = term_name.lower()
    matched = []
    for axis, keywords in axes_keywords.items():
        if any((kw in term_lower for kw in keywords)):
            matched.append(axis)
    return matched if matched else ['Other']
all_enrichment_del['axes'] = all_enrichment_del['Term'].apply(classify_axis)
enrichment_del_exploded = all_enrichment_del.explode('axes')
for axis_name in ['Cholinergic/Neurotransmitter', 'Neuroinflammatory', 'Neurodegenerative', 'Hormonal/Endocrine', 'Vascular/BBB']:
    subset = enrichment_del_exploded[enrichment_del_exploded['axes'] == axis_name]
    if len(subset) > 0:
        print(f"\n{'=' * 60}")
        print(f'Top {axis_name} terms:')
        print(f"{'=' * 60}")
        print(subset.nsmallest(5, 'Adjusted P-value')[['Term', 'Adjusted P-value', 'Genes', 'database']].to_string(index=False))
