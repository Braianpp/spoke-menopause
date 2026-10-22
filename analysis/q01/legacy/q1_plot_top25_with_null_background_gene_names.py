# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""Q1: Top-25 genes (Premature Menopause → Dementia) with null background overlay,
using *gene names* (from the SPOKE graph node attribute `name`) instead of node IDs.

Real bars:
  - total_occurrences from pm_dementia_dwpc_recurring_genes_v2.csv
Null background:
  - bg_mean ± bg_std from pm_dementia_null_gene_pvalues_v2.csv

Important:
- Null counts are non-negative integers, so mean-std can be negative.
  We plot asymmetric error bars with the *left* side clipped at 0.

Outputs an SVG with a q1_ prefix.
"""
from __future__ import annotations
import os
import argparse
import pickle
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
PM_COLOR = '#93C4BC'

def load_graph(gfile: Path):
    if not gfile.exists():
        raise FileNotFoundError(f'Graph not found: {gfile}')
    with open(gfile, 'rb') as fh:
        return pickle.load(fh)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--real-csv', default='${SPOKE_DATA_ROOT}/pm_dementia_dwpc_recurring_genes_v2.csv')
    ap.add_argument('--null-csv', default='${SPOKE_DATA_ROOT}/pm_dementia_null_gene_pvalues_v2.csv')
    ap.add_argument('--graph', default='${SPOKE_DATA_ROOT}/V1_homosapiens_1N_9E.gpickle', help='SPOKE graph gpickle used to map Gene:#### -> name')
    ap.add_argument('--top-n', type=int, default=25)
    ap.add_argument('--out', default='${SPOKE_DATA_ROOT}/biorouter/figs/q1_top25_genes_with_null_background_gene_names.svg')
    args = ap.parse_args()
    real_csv = Path(args.real_csv)
    null_csv = Path(args.null_csv)
    if not real_csv.exists():
        raise SystemExit(f'Missing real CSV: {real_csv}')
    if not null_csv.exists():
        raise SystemExit(f'Missing null CSV: {null_csv}')
    rec = pd.read_csv(real_csv)
    nul = pd.read_csv(null_csv)
    gcol = 'gene_name' if 'gene_name' in rec.columns else rec.columns[0]
    if 'total_occurrences' in rec.columns:
        real_col = 'total_occurrences'
    elif 'real_count' in rec.columns:
        real_col = 'real_count'
    else:
        raise SystemExit('Real CSV must have total_occurrences (or real_count)')
    for c in ['bg_mean', 'bg_std']:
        if c not in nul.columns:
            raise SystemExit(f'Null CSV missing required column: {c}')
    df = rec[[gcol, real_col]].merge(nul[[c for c in [gcol, 'bg_mean', 'bg_std', 'empirical_pval', 'significant_bh'] if c in nul.columns]], on=gcol, how='left')
    G = load_graph(Path(args.graph))

    def gid_to_name(gid: str) -> str:
        try:
            return G.nodes[gid].get('name', gid)
        except Exception:
            return gid
    df['gene_label'] = df[gcol].astype(str).map(gid_to_name)
    top = df.sort_values(real_col, ascending=False).head(args.top_n).copy()
    sns.set_theme(context='talk', style='whitegrid')
    plt.rcParams.update({'svg.fonttype': 'none'})
    fig, ax = plt.subplots(figsize=(8.0, 7.2))
    sns.barplot(data=top, y='gene_label', x=real_col, color=PM_COLOR, edgecolor='black', linewidth=0.6, ax=ax)
    y_pos = np.arange(len(top))
    bg_mean = pd.to_numeric(top['bg_mean'], errors='coerce').fillna(0.0).to_numpy()
    bg_std = pd.to_numeric(top['bg_std'], errors='coerce').fillna(0.0).to_numpy()
    xerr_left = np.minimum(bg_std, bg_mean)
    xerr_right = bg_std
    ax.errorbar(x=bg_mean, y=y_pos, xerr=[xerr_left, xerr_right], fmt='o', color='black', ecolor='black', elinewidth=1.5, capsize=3, label='Null mean ± SD (clipped at 0)', zorder=10)
    ax.set_title('Q1 — Top 25 Genes by Metapath Recurrence\n(Premature Menopause → Dementia)', fontsize=16, fontweight='bold')
    ax.set_xlabel('Real metapath-occurrence count', fontsize=13, fontweight='bold')
    ax.set_ylabel('')
    for lab in ax.get_yticklabels():
        lab.set_fontsize(11.5)
        lab.set_fontweight('bold')
    ax.legend(frameon=True, framealpha=0.95, edgecolor='#333333')
    sns.despine(ax=ax)
    fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, format='svg', bbox_inches='tight')
    print(f'Saved: {out}')
if __name__ == '__main__':
    main()
