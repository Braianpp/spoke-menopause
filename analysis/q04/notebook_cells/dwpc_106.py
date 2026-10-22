# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
import gseapy as gp
bone_brain_genes = list(bone_brain)
print(f'Running enrichment on {len(bone_brain_genes)} bone-brain overlap genes')
gene_set_libraries = ['KEGG_2021_Human', 'Reactome_2022', 'GO_Biological_Process_2023', 'WikiPathway_2023_Human', 'Jensen_DISEASES']
enrichment_results = {}
for lib in gene_set_libraries:
    try:
        enr = gp.enrichr(gene_list=bone_brain_genes, gene_sets=lib, organism='human', outdir=None, cutoff=0.05)
        enrichment_results[lib] = enr.results[enr.results['Adjusted P-value'] < 0.05]
        print(f'  {lib}: {len(enrichment_results[lib])} significant terms')
    except Exception as e:
        print(f'  {lib}: FAILED - {e}')
all_enrichment_bb = pd.concat([df.assign(database=db) for db, df in enrichment_results.items()], ignore_index=True)
all_enrichment_bb.to_csv('bone_brain_enrichment.csv', index=False)
print(f'\nTotal significant terms: {len(all_enrichment_bb)}')
axes_keywords = {'Estrogen/Hormonal': ['estrogen', 'steroid', 'hormone', 'gonad', 'ovari', 'androgen', 'aromatase', 'progesterone', 'menopause'], 'Bone/Skeletal': ['bone', 'osteo', 'osteoblast', 'osteoclast', 'skeletal', 'mineral', 'calcium', 'fracture', 'cartilage', 'BMP'], 'Inflammatory': ['inflam', 'cytokine', 'interleukin', 'TNF', 'NF-kB', 'chemokine', 'immune', 'toll-like', 'macrophage'], 'Neurodegenerative': ['neuro', 'synap', 'amyloid', 'tau', 'dopamin', 'axon', 'myelin', 'cognitive', 'brain', 'dementia', 'alzheimer'], 'Vascular': ['vascul', 'angiogen', 'endotheli', 'coagul', 'platelet', 'atheroscl', 'hypoxia']}

def classify_axis(term_name):
    term_lower = term_name.lower()
    matched = []
    for axis, keywords in axes_keywords.items():
        if any((kw in term_lower for kw in keywords)):
            matched.append(axis)
    return matched if matched else ['Other']
all_enrichment_bb['axes'] = all_enrichment_bb['Term'].apply(classify_axis)
enrichment_exploded = all_enrichment_bb.explode('axes')
for axis_name in ['Estrogen/Hormonal', 'Bone/Skeletal', 'Inflammatory', 'Neurodegenerative', 'Vascular']:
    subset = enrichment_exploded[enrichment_exploded['axes'] == axis_name]
    if len(subset) > 0:
        print(f"\n{'=' * 60}")
        print(f'Top {axis_name} terms:')
        print(f"{'=' * 60}")
        print(subset.nsmallest(5, 'Adjusted P-value')[['Term', 'Adjusted P-value', 'Genes', 'database']].to_string(index=False))
