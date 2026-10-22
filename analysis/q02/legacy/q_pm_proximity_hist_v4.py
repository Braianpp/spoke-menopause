# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""
Network Proximity: Premature Menopause → All Diseases in SPOKE (v4)

Reads pm_all_diseases_3hop.csv and produces a cleaner, publication-style
histogram with annotated reference lines for selected targets, optional
reference lines for other neurodegenerative diseases, and a separate
transparent legend export.

v4 changes
----------
* Alzheimer's disease is promoted to its own highlighted line (distinct
  colour, own legend entry) instead of being one of the anonymous purple
  "other neurodegenerative" lines.
* Reference lines whose scores collide (identical or within a tolerance)
  are drawn as side-by-side horizontal segments rather than stacked on
  top of one another, so no line can be completely hidden.
* Optional --annotate-neuro to label the highlighted/other neuro lines
  directly on the plot.

Usage:
  python q_pm_proximity_hist_v4.py       --csv pm_all_diseases_3hop.csv       --outdir biorouter/figs       --dpi 600       --bin-size auto       --include-other-neuro
"""
import os
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'axes.linewidth': 0.8})
C_HIST = '#D9D9D9'
C_TEXT = '#333333'
C_SUB = '#666666'
C_OTHER_NEURO = '#9B59B6'
DISEASE_COLORS = {'Vascular dementia': '#E74C3C', 'Frontotemporal dementia': '#3498DB', "Parkinson's disease": '#2ECC71', "Alzheimer's disease": '#E67E22'}
SELECTED_TARGETS = {'DOID:8725': 'Vascular dementia', 'DOID:9255': 'Frontotemporal dementia', 'DOID:14330': "Parkinson's disease"}
HIGHLIGHT_NEURO = {'DOID:10652': "Alzheimer's disease"}
OTHER_NEURO = {'DOID:332': 'Amyotrophic lateral sclerosis', 'DOID:12858': "Huntington's disease", 'DOID:678': 'Multiple sclerosis', 'DOID:1289': 'Neurodegenerative disease', 'DOID:12217': 'Lewy body dementia', 'DOID:2474': 'Dementia'}

def ensure_outdir(p: str):
    os.makedirs(p, exist_ok=True)

def auto_bin_size(values: np.ndarray) -> int:
    vmax = int(np.nanmax(values))
    if vmax <= 40:
        return 1
    if vmax <= 100:
        return 2
    return 5

def save_transparent_legend(handles, labels, outpath, ncol=1, fontsize=10, title=None):
    fig = plt.figure(figsize=(2.5, 2.0))
    fig.patch.set_alpha(0.0)
    ax = fig.add_subplot(111)
    ax.axis('off')
    fig.legend(handles=handles, labels=labels, loc='center', frameon=False, ncol=ncol, fontsize=fontsize, title=title, title_fontsize=fontsize)
    fig.savefig(outpath, dpi=600, bbox_inches='tight', transparent=True)
    plt.close(fig)

def load_df(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    need = {'disease_id', 'disease', 'shared_3hop'}
    if not need.issubset(set(df.columns)):
        raise ValueError(f'Missing required columns in {csv_path}: {need - set(df.columns)}')
    df = df.sort_values('shared_3hop', ascending=False).reset_index(drop=True)
    if 'rank' not in df.columns:
        df['rank'] = np.arange(1, len(df) + 1)
    return df

def _score_for(df: pd.DataFrame, disease_id: str):
    row = df[df['disease_id'] == disease_id]
    if row.empty:
        return None
    return float(row.iloc[0]['shared_3hop'])

def draw_level_lines(ax, entries, tol_frac=0.005, annotate=False, gap=0.012):
    """Draw horizontal reference lines, splitting collisions side by side.

    entries: list of dicts with keys score, label, color, lw, ls, zorder, alpha.
    Lines whose scores fall within tol of each other are drawn as adjacent
    horizontal segments spanning different fractions of the axes width, so
    every line stays visible even when the underlying values are identical.
    """
    if not entries:
        return
    ymin, ymax = ax.get_ylim()
    tol = tol_frac * (ymax - ymin)
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    order = sorted(range(len(entries)), key=lambda i: entries[i]['score'])
    groups, cur = ([], [order[0]])
    for i in order[1:]:
        if abs(entries[i]['score'] - entries[cur[-1]]['score']) <= tol:
            cur.append(i)
        else:
            groups.append(cur)
            cur = [i]
    groups.append(cur)
    for g in groups:
        n = len(g)
        g = sorted(g, key=lambda i: -entries[i]['zorder'])
        for k, idx in enumerate(g):
            e = entries[idx]
            if n == 1:
                x0, x1 = (0.0, 1.0)
            else:
                x0 = k / n
                x1 = (k + 1) / n - gap
            ax.plot([x0, x1], [e['score'], e['score']], transform=trans, color=e['color'], linewidth=e['lw'], linestyle=e['ls'], alpha=e.get('alpha', 1.0), zorder=e['zorder'], solid_capstyle='butt')
            if annotate and e.get('annotate', False):
                ax.text(min(x1, 0.995), e['score'] + 0.006 * (ymax - ymin), e['label'], transform=trans, ha='right', va='bottom', fontsize=8, color=e['color'], zorder=e['zorder'] + 1)

def plot_proximity(df: pd.DataFrame, out_png: str, out_pdf: str, out_legend: str | None=None, bin_size: str | int='auto', include_other_neuro: bool=False, dpi: int=600, xmax_override: float | None=None, include_callouts: bool=False, show_legend_in_figure: bool=False, show_stats: bool=False, annotate_neuro: bool=False):
    positive = df[df['shared_3hop'] > 0].copy()
    if positive.empty:
        raise ValueError('No diseases with shared_3hop > 0.')
    vmax = int(positive['shared_3hop'].max())
    bin_w = auto_bin_size(positive['shared_3hop']) if bin_size == 'auto' else int(bin_size)
    bins = np.arange(0, vmax + bin_w + 1, bin_w)
    fig, ax = plt.subplots(figsize=(6.0, 8.0))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, color='#ECECEC', linewidth=0.7)
    ax.hist(positive['shared_3hop'], bins=bins, color=C_HIST, edgecolor='#FFFFFF', linewidth=0.6, alpha=0.95, zorder=2, orientation='horizontal')
    ymax = vmax + max(3, int(0.05 * vmax))
    if xmax_override is not None:
        ymax = float(xmax_override)
    ax.set_ylim(0, ymax)
    ax.set_ylabel('Shared 3‑Hop Gene Intermediates (PM → Gene → Gene → Disease)', fontsize=12, fontweight='bold', labelpad=6)
    ax.set_xlabel('Number of Diseases', fontsize=12, fontweight='bold', labelpad=6)
    ax.tick_params(axis='both', labelsize=11)
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    ax.set_title('Network Proximity: Premature Menopause → All Diseases in SPOKE', fontsize=14, fontweight='bold', loc='left', pad=10)
    entries = []
    legend_lines = []
    for did, dname in SELECTED_TARGETS.items():
        score = _score_for(df, did)
        if score is None:
            continue
        color = DISEASE_COLORS[dname]
        entries.append(dict(score=score, label=dname, color=color, lw=3.0, ls='--', zorder=6, alpha=0.95, annotate=False))
        legend_lines.append(Line2D([0], [0], color=color, linewidth=3.0, linestyle='--', label=dname))
    for did, dname in HIGHLIGHT_NEURO.items():
        score = _score_for(df, did)
        if score is None or score <= 0:
            continue
        color = DISEASE_COLORS.get(dname, C_OTHER_NEURO)
        entries.append(dict(score=score, label=dname, color=color, lw=2.6, ls=(0, (5, 2)), zorder=5, alpha=1.0, annotate=True))
        legend_lines.append(Line2D([0], [0], color=color, linewidth=2.6, linestyle=(0, (5, 2)), label=dname))
    if include_other_neuro:
        for did, name in OTHER_NEURO.items():
            score = _score_for(df, did)
            if score is None or score <= 0:
                continue
            entries.append(dict(score=score, label=name, color=C_OTHER_NEURO, lw=1.8, ls='-', zorder=3, alpha=0.65, annotate=True))
    draw_level_lines(ax, entries, annotate=annotate_neuro)
    if include_callouts:
        ty = blended_transform_factory(ax.transAxes, ax.transData)
        for i, (did, dname) in enumerate(SELECTED_TARGETS.items()):
            score = _score_for(df, did)
            if score is None:
                continue
            color = DISEASE_COLORS[dname]
            x_off = [1.05, 1.13, 1.21][i % 3]
            ax.annotate(dname, xy=(1.0, score), xycoords=ty, xytext=(x_off, score), textcoords=ty, va='center', ha='left', fontsize=10, fontweight='bold', color=color, bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor=color, linewidth=1.0, alpha=0.98), arrowprops=dict(arrowstyle='-', color=color, linewidth=1.0), annotation_clip=False, zorder=7)
    legend_handles = [Patch(facecolor=C_HIST, edgecolor='#FFFFFF', label='All diseases')] + legend_lines
    if include_other_neuro:
        legend_handles.append(Line2D([0], [0], color=C_OTHER_NEURO, linewidth=1.8, linestyle='-', alpha=0.65, label='Other neurodegenerative'))
    legend_labels = [h.get_label() for h in legend_handles]
    if show_legend_in_figure:
        ax.legend(handles=legend_handles, labels=legend_labels, loc='upper right', fontsize=9, frameon=True, framealpha=0.96, edgecolor='#DDDDDD', fancybox=True)
    if show_stats:
        yvals = positive['shared_3hop'].values
        stats_text = f'Diseases: {len(df):,}\nConnected: {(yvals > 0).sum()}\nMedian: {np.median(yvals):.0f}\nMean: {np.mean(yvals):.1f}'
        ax.text(0.99, 0.03, stats_text, transform=ax.transAxes, ha='right', va='bottom', fontsize=9, color=C_TEXT, bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='#DDDDDD', alpha=0.96))
    plt.tight_layout()
    fig.savefig(out_png, dpi=dpi, bbox_inches='tight', facecolor='white')
    fig.savefig(out_pdf, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    if out_legend:
        save_transparent_legend(legend_handles, legend_labels, out_legend, ncol=1, fontsize=10, title=None)

def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--csv', required=True, help='pm_all_diseases_3hop.csv')
    ap.add_argument('--outdir', default='.', help='Output directory')
    ap.add_argument('--dpi', type=int, default=600, help='DPI for saved figures')
    ap.add_argument('--bin-size', default='auto', help="Histogram bin size (int or 'auto')")
    ap.add_argument('--include-other-neuro', action='store_true', help='Show de-emphasized lines for other neurodegenerative diseases')
    ap.add_argument('--annotate-neuro', action='store_true', help='Label the neurodegenerative reference lines on the plot')
    ap.add_argument('--xmax', type=float, default=180.0, help='Maximum y-axis (shared intermediates)')
    args = ap.parse_args()
    ensure_outdir(args.outdir)
    df = load_df(args.csv)
    out_png = os.path.join(args.outdir, 'pm_all_diseases_3hop_proximity_v4.png')
    out_pdf = os.path.join(args.outdir, 'pm_all_diseases_3hop_proximity_v4.pdf')
    out_legend = os.path.join(args.outdir, 'pm_all_diseases_3hop_legend_transparent.png')
    plot_proximity(df, out_png, out_pdf, out_legend=out_legend, bin_size=args.bin_size, include_other_neuro=args.include_other_neuro, dpi=args.dpi, xmax_override=args.xmax, annotate_neuro=args.annotate_neuro)
    print('Saved:\n ', out_png, '\n ', out_pdf, '\n ', out_legend)
if __name__ == '__main__':
    main()
