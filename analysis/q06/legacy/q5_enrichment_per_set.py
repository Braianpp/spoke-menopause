# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import gseapy as gp
GENE_CSV = '${SPOKE_DATA_ROOT}/dwpc_delirium_genes.csv'
OUT_DIR = '${SPOKE_DATA_ROOT}/biorouter/figs'
CSV_DIR = '${SPOKE_DATA_ROOT}/biorouter'
PAIR_PM_DEL = 'Premature menopause → Subacute delirium'
PAIR_PM_DEM = 'Premature menopause → Dementia'
PAIR_DEL_DEM = 'Subacute delirium → Dementia'
PAIR_PM_DEL_ASCII = 'Premature menopause -> Subacute delirium'
PAIR_PM_DEM_ASCII = 'Premature menopause -> Dementia'
PAIR_DEL_DEM_ASCII = 'Subacute delirium -> Dementia'
libraries = ['KEGG_2021_Human', 'Reactome_2022', 'GO_Biological_Process_2023', 'WikiPathway_2023_Human', 'Jensen_DISEASES']
axis_colors = {'Cholinergic/Neurotransmitter': '#F39C12', 'Neuroinflammatory': '#E67E22', 'Neurodegenerative': '#2C3E50', 'Hormonal/Endocrine': '#9B59B6', 'Vascular/BBB': '#C0392B', 'Other': '#BDC3C7'}
axes_keywords = {'Cholinergic/Neurotransmitter': ['acetylchol', 'cholinerg', 'dopamin', 'serotonin', 'gaba', 'glutamat', 'neurotransmit', 'synapt'], 'Neuroinflammatory': ['inflam', 'cytokine', 'interleukin', 'tnf', 'nf-kb', 'chemokine', 'immune', 'toll-like', 'microglia', 'astrocyte'], 'Neurodegenerative': ['neuro', 'amyloid', 'tau', 'dementia', 'alzheimer', 'cognitive', 'axon', 'myelin'], 'Hormonal/Endocrine': ['estrogen', 'steroid', 'hormone', 'gonad', 'ovari', 'androgen', 'aromatase', 'menopause'], 'Vascular/BBB': ['vascul', 'angiogen', 'endotheli', 'blood-brain', 'hypoxia', 'ischem', 'stroke', 'cerebrovascul']}

def classify_axis(term: str) -> str:
    t = str(term).lower()
    for ax_name, kws in axes_keywords.items():
        if any((kw in t for kw in kws)):
            return ax_name
    return 'Other'

def get_genes(df: pd.DataFrame, label: str) -> list:
    mask = df['pairs'].fillna('').astype(str).str.contains(label)
    if not mask.any():
        mask = df['pairs'].fillna('').astype(str).str.contains(label.replace('→', '->'))
    genes = sorted(set(df.loc[mask, 'gene_name'].dropna().astype(str)))
    return genes

def main():
    if not os.path.exists(GENE_CSV):
        raise SystemExit(f'Gene CSV not found on leo: {GENE_CSV}')
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(CSV_DIR, exist_ok=True)
    df = pd.read_csv(GENE_CSV)
    if 'pair' in df.columns and 'pairs' not in df.columns:
        df = df.rename(columns={'pair': 'pairs'})
    for required in ['pairs', 'gene_name']:
        if required not in df.columns:
            raise SystemExit(f"Missing column '{required}' in {GENE_CSV}")
    g_pm_del = get_genes(df, PAIR_PM_DEL)
    g_pm_dem = get_genes(df, PAIR_PM_DEM)
    g_del_dem = get_genes(df, PAIR_DEL_DEM)
    print(f'PM→Delirium genes: {len(g_pm_del)}')
    print(f'PM→Dementia genes: {len(g_pm_dem)}')
    print(f'Delirium→Dementia genes: {len(g_del_dem)}')
    results = {}
    for label, genes in [(PAIR_PM_DEM, g_pm_dem), (PAIR_DEL_DEM, g_del_dem), (PAIR_PM_DEL, g_pm_del)]:
        if len(genes) < 3:
            print(f'Skipping enrichment for {label} (n={len(genes)})')
            results[label] = pd.DataFrame(columns=['Term', 'Adjusted P-value', 'Genes', 'database', 'pair'])
            continue
        all_rows = []
        for lib in libraries:
            try:
                enr = gp.enrichr(gene_list=genes, gene_sets=lib, organism='human', outdir=None, cutoff=0.05)
                df_en = enr.results
                if df_en is not None and (not df_en.empty):
                    df_en = df_en[df_en['Adjusted P-value'] < 0.05].copy()
                    if not df_en.empty:
                        df_en['database'] = lib
                        df_en['pair'] = label
                        all_rows.append(df_en)
                        print(f'  {label} | {lib}: {len(df_en)} sig terms')
            except Exception as e:
                print(f'  {label} | {lib}: FAILED - {e}')
        if all_rows:
            res = pd.concat(all_rows, ignore_index=True)
        else:
            res = pd.DataFrame(columns=['Term', 'Adjusted P-value', 'Genes', 'database', 'pair'])
        results[label] = res
        safe = label.replace(' → ', '_to_').replace(' ', '_')
        res.to_csv(os.path.join(CSV_DIR, f'q5_enrichment_{safe}.csv'), index=False)
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle('Q5: Enrichment per gene set (DWPC-derived)', fontsize=14, fontweight='bold', y=1.02)
    for ax, label in zip(axes, [PAIR_PM_DEM, PAIR_DEL_DEM]):
        res = results.get(label, pd.DataFrame())
        ax.set_title(label, fontsize=12, fontweight='bold')
        if res is None or res.empty:
            ax.text(0.5, 0.5, 'No significant terms', ha='center', va='center')
            continue
        res = res.copy()
        res['axis'] = res['Term'].map(classify_axis)
        cur = res.nsmallest(10, 'Adjusted P-value')
        names = cur['Term'].str.slice(0, 60).tolist()[::-1]
        vals = (-np.log10(cur['Adjusted P-value']).values)[::-1]
        cols = [axis_colors.get(a, '#BDC3C7') for a in cur['axis'][::-1]]
        ax.barh(names, vals, color=cols, edgecolor='white', linewidth=0.5, height=0.6)
        ax.set_xlabel('-log10(Adjusted P-value)')
        ax.axvline(x=-np.log10(0.05), color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    legend_patches = [mpatches.Patch(color=c, label=a) for a, c in axis_colors.items() if a != 'Other']
    axes[1].legend(handles=legend_patches, loc='lower right', fontsize=8)
    plt.tight_layout()
    figpath = os.path.join(OUT_DIR, 'q5_enrichment_per_set.png')
    plt.savefig(figpath, dpi=300, bbox_inches='tight')
    print('Saved:', figpath)
if __name__ == '__main__':
    main()
