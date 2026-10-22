# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
import os
import time
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
LIBS_PATHWAY = ['KEGG_2021_Human', 'Reactome_2022', 'GO_Biological_Process_2023', 'WikiPathway_2023_Human']
AXIS_COLORS = {'Cholinergic/Neurotransmitter': '#F39C12', 'Neuroinflammatory': '#E67E22', 'Neurodegenerative': '#2C3E50', 'Hormonal/Endocrine': '#9B59B6', 'Vascular/BBB': '#C0392B', 'Other': '#BDC3C7'}
AXIS_KEYWORDS = {'Cholinergic/Neurotransmitter': ['acetylchol', 'cholinerg', 'dopamin', 'serotonin', 'gaba', 'glutamat', 'neurotransmit', 'synapt'], 'Neuroinflammatory': ['inflam', 'cytokine', 'interleukin', 'tnf', 'nf-kb', 'chemokine', 'immune', 'toll-like', 'microglia', 'astrocyte'], 'Neurodegenerative': ['neuro', 'amyloid', 'tau', 'alzheimer', 'cognitive', 'axon', 'myelin'], 'Hormonal/Endocrine': ['estrogen', 'steroid', 'hormone', 'gonad', 'ovari', 'androgen', 'aromatase', 'menopause'], 'Vascular/BBB': ['vascul', 'angiogen', 'endotheli', 'blood-brain', 'hypoxia', 'ischem', 'stroke', 'cerebrovascul']}
EXCLUDE_TERMS = ['disease', 'dementia', 'encephalopathy']

def classify_axis(term: str) -> str:
    t = str(term).lower()
    for ax_name, kws in AXIS_KEYWORDS.items():
        if any((kw in t for kw in kws)):
            return ax_name
    return 'Other'

def get_genes(df: pd.DataFrame, label: str) -> list:
    s = df['pairs'].fillna('').astype(str)
    mask = s.str.contains(label) | s.str.contains(label.replace('→', '->'))
    genes = sorted(set(df.loc[mask, 'gene_name'].dropna().astype(str)))
    return genes

def enrich_with_retry(genes: list, gene_set: str, max_retries: int=4, base_wait: float=1.0):
    for attempt in range(max_retries):
        try:
            enr = gp.enrichr(gene_list=genes, gene_sets=gene_set, organism='human', outdir=None, cutoff=0.05)
            return enr
        except Exception as e:
            if attempt == max_retries - 1:
                raise
            time.sleep(base_wait * 2 ** attempt)
    return None

def main():
    if not os.path.exists(GENE_CSV):
        raise SystemExit(f'Gene CSV not found: {GENE_CSV}')
    df = pd.read_csv(GENE_CSV)
    if 'pair' in df.columns and 'pairs' not in df.columns:
        df = df.rename(columns={'pair': 'pairs'})
    for col in ['pairs', 'gene_name']:
        if col not in df.columns:
            raise SystemExit(f"Missing column '{col}' in {GENE_CSV}")
    g_pm_del = get_genes(df, PAIR_PM_DEL)
    g_pm_dem = get_genes(df, PAIR_PM_DEM)
    g_del_dem = get_genes(df, PAIR_DEL_DEM)
    print(f'PM→Delirium genes: {len(g_pm_del)}')
    print(f'PM→Dementia genes: {len(g_pm_dem)}')
    print(f'Delirium→Dementia genes: {len(g_del_dem)}')
    pair_lists = [(PAIR_PM_DEM, g_pm_dem), (PAIR_DEL_DEM, g_del_dem), (PAIR_PM_DEL, g_pm_del)]
    results = {}
    for label, genes in pair_lists:
        if len(genes) < 3:
            print(f'Skipping {label} (n={len(genes)})')
            results[label] = pd.DataFrame(columns=['Term', 'Adjusted P-value', 'Genes', 'database', 'pair'])
            continue
        rows = []
        for lib in LIBS_PATHWAY:
            try:
                enr = enrich_with_retry(genes, lib)
                df_en = enr.results if enr is not None else pd.DataFrame()
                if df_en is not None and (not df_en.empty):
                    df_en = df_en[df_en['Adjusted P-value'] < 0.05].copy()
                    df_en = df_en[~df_en['Term'].str.lower().str.contains('|'.join(EXCLUDE_TERMS))]
                    if not df_en.empty:
                        df_en['database'] = lib
                        df_en['pair'] = label
                        rows.append(df_en)
                        print(f'  {label} | {lib}: {len(df_en)} sig pathways')
            except Exception as e:
                print(f'  {label} | {lib}: FAILED - {e}')
        res = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=['Term', 'Adjusted P-value', 'Genes', 'database', 'pair'])
        results[label] = res
        safe = label.replace(' → ', '_to_').replace(' ', '_')
        res.to_csv(os.path.join(CSV_DIR, f'q5_pathway_enrichment_{safe}.csv'), index=False)
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle('Q5: Pathway enrichment per gene set (DWPC-derived)', fontsize=14, fontweight='bold', y=1.02)
    for ax, label in zip(axes, [PAIR_PM_DEM, PAIR_PM_DEL]):
        res = results.get(label, pd.DataFrame())
        ax.set_title(label, fontsize=12, fontweight='bold')
        if res is None or res.empty:
            ax.text(0.5, 0.5, 'No significant pathways', ha='center', va='center')
            continue
        res = res.copy()
        res['axis'] = res['Term'].map(classify_axis)
        cur = res.nsmallest(10, 'Adjusted P-value')
        names = cur['Term'].str.slice(0, 60).tolist()[::-1]
        vals = (-np.log10(cur['Adjusted P-value']).values)[::-1]
        cols = [AXIS_COLORS.get(a, '#BDC3C7') for a in cur['axis'][::-1]]
        ax.barh(names, vals, color=cols, edgecolor='white', linewidth=0.5, height=0.6)
        ax.set_xlabel('-log10(Adjusted P-value)')
        ax.axvline(x=-np.log10(0.05), color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    legend_patches = [mpatches.Patch(color=c, label=a) for a, c in AXIS_COLORS.items() if a != 'Other']
    axes[1].legend(handles=legend_patches, loc='lower right', fontsize=8)
    plt.tight_layout()
    figpath = os.path.join(OUT_DIR, 'q5_pathway_enrichment_per_set.png')
    plt.savefig(figpath, dpi=300, bbox_inches='tight')
    print('Saved:', figpath)
if __name__ == '__main__':
    main()
