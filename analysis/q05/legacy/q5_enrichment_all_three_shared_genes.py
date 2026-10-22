# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""
Enrichment analysis (GO_Biological_Process_2023) on the 36 genes shared
by all three PPR embeddings (PM, APO, Cognitive) at top 2%.

Output:
  - CSV: q5_all_three_enrichr_GO_BP_2023.csv
  - SVG: q5_all_three_enrichment_GO_BP_2023.svg
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gseapy as gp
import os
GENES_CSV = '${SPOKE_DATA_ROOT}/shared_all_three_v4b_top2pct_2hops.csv'
OUT_DIR = '${SPOKE_DATA_ROOT}/biorouter/figs'
RESULTS_DIR = '${SPOKE_DATA_ROOT}/biorouter/results'
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
BAR_COLOR = '#93C4BC'
print('Loading shared gene list...')
df = pd.read_csv(GENES_CSV)
gene_symbols = sorted(df['gene_name'].dropna().astype(str).tolist())
print(f'Genes: {len(gene_symbols)}')
print(f"  {', '.join(gene_symbols[:10])}...")
print('\nRunning GO Biological Process 2023 enrichment...')
enr = gp.enrichr(gene_list=gene_symbols, gene_sets=['GO_Biological_Process_2023'], organism='human', outdir=None, no_plot=True)
enr_df = enr.results
enr_df_sig = enr_df[enr_df['Adjusted P-value'] < 0.05].copy()
print(f'Total enrichment terms: {len(enr_df)}')
print(f'Significant terms (FDR < 0.05): {len(enr_df_sig)}')
enr_csv_path = os.path.join(RESULTS_DIR, 'q5_all_three_enrichr_GO_BP_2023.csv')
enr_df.to_csv(enr_csv_path, index=False)
print(f'Saved enrichment results: {enr_csv_path}')
plot_df = enr_df.sort_values('Adjusted P-value').head(12).copy()
plot_df['-log10(FDR)'] = -np.log10(plot_df['Adjusted P-value'])
plot_df['Term_short'] = plot_df['Term'].apply(lambda x: x.split(' (')[0] if ' (' in x else x)
fig, ax = plt.subplots(figsize=(10, 7))
bars = ax.barh(range(len(plot_df)), plot_df['-log10(FDR)'].values, color=BAR_COLOR, edgecolor='black', linewidth=0.8)
ax.set_yticks(range(len(plot_df)))
ax.set_yticklabels(plot_df['Term_short'].values, fontsize=12, fontweight='bold')
ax.invert_yaxis()
ax.set_xlabel('-log₁₀(Adjusted P-value)', fontsize=14, fontweight='bold')
ax.set_title('GO Biological Process 2023 Enrichment\nShared Genes (PM ∩ APO ∩ Cognitive, top 2%)', fontsize=14, fontweight='bold')
ax.axvline(-np.log10(0.05), color='red', linestyle='--', linewidth=1, alpha=0.7, label='FDR = 0.05')
ax.legend(loc='lower right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.yaxis.grid(False)
ax.xaxis.grid(True, color='gray', alpha=0.3, linestyle='-', linewidth=0.5)
plt.tight_layout()
enr_fig_path = os.path.join(OUT_DIR, 'q5_all_three_enrichment_GO_BP_2023.svg')
fig.savefig(enr_fig_path, format='svg', bbox_inches='tight', dpi=150)
plt.close()
print(f'Saved enrichment plot: {enr_fig_path}')
print('\nDone!')
