# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
from pathlib import Path
import json, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
P = Path(__file__).resolve().parent
ROOT = P.parent
OUT = ROOT / 'figures/figure8_2026-10-01'
OUT.mkdir(exist_ok=True)
p = json.loads((ROOT / 'supplementary_tables_2026-10-01/payload_v3.json').read_text())
s = next((s for s in p['sections'] if s['name'] == 'S9g_Letrozole_metapaths'))
data = [dict(zip(s['headers'], r)) for r in s['rows'] if r[-1]]
assert len(data) == 10
plt.rcParams.update({'font.family': 'Arial', 'font.size': 17, 'text.color': '#233E4D', 'axes.labelcolor': '#233E4D', 'xtick.color': '#233E4D', 'ytick.color': '#233E4D', 'svg.fonttype': 'none', 'pdf.fonttype': 42})
fig, axes = plt.subplots(2, 1, figsize=(15, 10), gridspec_kw={'hspace': 0.54})
for ax, length in zip(axes, [3, 4]):
    rows = sorted([r for r in data if r['length'] == length], key=lambda r: -r['dwpc'])
    labs = [r['metapath'].replace(' → ', ' →\n', 1) if length == 4 else r['metapath'] for r in rows]
    labs = [' → '.join((x.split('_')[-1] for x in r['metapath'].split(' → '))) for r in rows]
    ax.barh(range(5), [r['dwpc'] for r in rows], height=0.6, color='#99BDB8', edgecolor='#698880', linewidth=0.8)
    ax.set_yticks(range(5), labs, fontsize=18)
    ax.invert_yaxis()
    ax.set_xlim(0, 0.00365)
    ax.set_title('A. Three-edge metapaths' if length == 3 else 'B. Four-edge metapaths', loc='left', color='black', fontsize=23, pad=16)
    ax.set_xlabel('DWPC score', fontsize=19, labelpad=8)
    sf = ScalarFormatter(useMathText=True)
    sf.set_powerlimits((-3, -3))
    ax.xaxis.set_major_formatter(sf)
    ax.grid(axis='x', color='#E5E9EC', linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.spines['bottom'].set_color('#95A6B0')
    ax.tick_params(axis='y', length=0, pad=10)
    for i, r in enumerate(rows):
        ax.text(r['dwpc'] + 4.5e-05, i, f"{r['dwpc'] * 1000:.3f}", va='center', fontsize=16)
fig.subplots_adjust(left=0.39, right=0.96, top=0.92, bottom=0.2)
fig.text(0.035, 0.085, 'CdG  Compound downregulates gene     CuG  Compound upregulates gene     DaG  Disease–gene association\nGPdG / KGdG  Gene downregulation     GPuG / KGuG  Gene upregulation\nCiC  Compound–compound interaction     PiC  Protein–compound interaction     PrG  Protein regulates gene', fontsize=13, ha='left', va='center', linespacing=1.5)
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(OUT / f'Figure8_letrozole_metapaths_corrected.{ext}', dpi=300, facecolor='white')
fig.savefig(OUT / 'Figure8_letrozole_metapaths_preview.png', dpi=90, facecolor='white')
plt.close(fig)
with (OUT / 'Figure8_plotted_metapaths.csv').open('w') as f:
    w = csv.DictWriter(f, fieldnames=list(data[0]))
    w.writeheader()
    w.writerows(data)
(OUT / 'Figure8_metapaths_caption.txt').write_text('Letrozole–dementia metapaths in the archived V1 graph. The five largest DWPC scores are shown separately for schemas with three edges (a) and four edges (b). Bars and labels use a common scale; numeric annotations are DWPC × 10³. Predicate codes identify relationship types and their original SPOKE prefixes. This is a redraw of saved values, not a new DWPC computation. Longer schemas permit different walks, so their scores do not establish a stronger biological effect. Source: Supplementary Table 9g.\n')
print(OUT)
