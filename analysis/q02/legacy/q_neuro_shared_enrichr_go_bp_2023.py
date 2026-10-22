# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Enrichr GO Biological Process enrichment for genes shared by >=2 neurodegenerative targets.

Gene set: genes with n_targets >= 2 in
  ${SPOKE_DATA_ROOT}/dwpc_recurring_genes.csv
(columns: node_id, gene_name, n_targets, n_metapaths, total_occurrences, targets, max_dwpc)

Method mirrors the Q1/Q3-era pipeline:
- Enrichr library GO_Biological_Process_2023 via gseapy
- BH-adjusted p-values
- Top-20 barplot (-log10 adj p), n= overlap labels

Outputs (q_neuro_shared prefix):
- CSV: ${SPOKE_DATA_ROOT}/biorouter/results/q_neuro_shared_enrichr_GO_Biological_Process_2023.csv
- SVG: ${SPOKE_DATA_ROOT}/biorouter/figs/q_neuro_shared_enrichr_GO_Biological_Process_2023_barplot.svg
"""
from __future__ import annotations
import os
import math
import re
from pathlib import Path
import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
LIBRARY = 'GO_Biological_Process_2023'
PM_COLOR = '#93C4BC'
GENES_CSV = Path('${SPOKE_DATA_ROOT}/dwpc_recurring_genes.csv')
OUT_RESULTS = Path(f'${{SPOKE_DATA_ROOT}}/biorouter/results/q_neuro_shared_enrichr_{LIBRARY}.csv')
OUT_PLOT = Path(f'${{SPOKE_DATA_ROOT}}/biorouter/figs/q_neuro_shared_enrichr_{LIBRARY}_barplot.svg')
TOP = 20
MIN_GENES = 2

def pick_adj_p_col(res):
    for c in ['Adjusted P-value', 'Adjusted P-value (Benjamini)', 'Adjusted P-value (BH)']:
        if c in res.columns:
            return c
    raise ValueError(f'No adj p col. Columns: {list(res.columns)}')

def parse_gene_list(s):
    if s is None:
        return []
    s = str(s)
    if s.lower() in {'nan', 'none'}:
        return []
    return [p.strip() for p in s.replace(';', ',').split(',') if p.strip()]

def main():
    df = pd.read_csv(GENES_CSV)
    df['n_targets'] = pd.to_numeric(df['n_targets'], errors='coerce')
    shared = df[df['n_targets'] >= 2].copy()
    symbols = sorted(set((str(g).strip() for g in shared['gene_name'] if str(g).strip() and str(g).lower() != 'nan')))
    print(f'[neuro-shared Enrichr] genes with n_targets>=2: {len(shared)}')
    print(f'[neuro-shared Enrichr] unique symbols: {len(symbols)}')
    print('Genes:', ', '.join(symbols))
    if len(symbols) < 5:
        raise SystemExit('Too few gene symbols for reliable enrichment.')
    enr = gp.enrichr(gene_list=symbols, gene_sets=[LIBRARY], organism='human', outdir=None, no_plot=True)
    res = enr.results
    if res is None or res.empty:
        raise SystemExit('No enrichment results returned by Enrichr.')
    res = res.copy()
    adj = pick_adj_p_col(res)
    if 'Genes' in res.columns:
        res['overlap_genes'] = res['Genes'].astype(str)
        res['n_overlap'] = res['overlap_genes'].map(lambda s: len(parse_gene_list(s)))
    else:
        res['overlap_genes'] = ''
        res['n_overlap'] = np.nan
    res_f = res[res['n_overlap'] >= MIN_GENES].copy()
    res_f = res_f.sort_values(adj, ascending=True).reset_index(drop=True)
    OUT_RESULTS.parent.mkdir(parents=True, exist_ok=True)
    res_f.to_csv(OUT_RESULTS, index=False)
    print(f'Saved results CSV: {OUT_RESULTS} (n={len(res_f)})')
    topn = min(TOP, len(res_f))
    plot_df = res_f.head(topn).copy()
    ABBREV = {'Positive Regulation Of Cysteine-Type Endopeptidase Activity Involved In Apoptotic Process': 'Pos. reg. apoptotic cysteine endopeptidase activity', 'Activation Of Cysteine-Type Endopeptidase Activity Involved In Apoptotic Process': 'Activation of apoptotic cysteine endopeptidase activity', 'Positive Regulation Of Protein Modification Process': 'Pos. reg. protein modification process'}

    def _clean_term(t):
        base = re.sub('\\s*\\(GO:\\d+\\)\\s*$', '', str(t)).strip()
        return ABBREV.get(base, base)
    plot_df['Term_label'] = plot_df['Term'].map(_clean_term)
    plot_df['minus_log10_adj_p'] = plot_df[adj].astype(float).map(lambda p: -math.log10(p) if p and p > 0 else np.nan)
    sns.set_theme(context='talk', style='whitegrid')
    plt.rcParams.update({'svg.fonttype': 'none'})
    fig_h = max(4.0, 0.35 * topn + 1.5)
    fig, ax = plt.subplots(figsize=(10.0, fig_h))
    sns.barplot(data=plot_df, y='Term_label', x='minus_log10_adj_p', color=PM_COLOR, edgecolor='black', linewidth=0.4, ax=ax)
    ax.set_title(f'Enrichr: {LIBRARY}\n(Genes shared by >=2 neurodegenerative targets: FTD/VaD/PD)', fontsize=14, fontweight='bold')
    ax.set_xlabel('-log10(adjusted p-value)', fontsize=12, fontweight='bold')
    ax.set_ylabel('')
    for i, r in plot_df.reset_index(drop=True).iterrows():
        try:
            ax.text(float(r['minus_log10_adj_p']) + 0.02, i, f"n={int(r['n_overlap'])}", va='center', fontsize=9)
        except Exception:
            pass
    sns.despine(ax=ax)
    fig.tight_layout()
    OUT_PLOT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PLOT, format='svg', bbox_inches='tight')
    fig.savefig(str(OUT_PLOT).replace('.svg', '.png'), dpi=300, bbox_inches='tight')
    print(f'Saved barplot: {OUT_PLOT} (+ .png)')
if __name__ == '__main__':
    main()
