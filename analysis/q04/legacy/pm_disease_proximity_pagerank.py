# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""
Generate the figure(s):
  "Disease Network Proximity to Premature Menopause (PageRank embedding from PM node)"

This standalone script:
- Loads a personalized PageRank embedding seeded at the Premature Menopause (PM) disease node
- Filters to Disease nodes, ranks by proximity (PageRank), computes percentiles
- Classifies diseases into coarse categories for visualization
- Outputs:
  1) Histogram panel only (legend removed), with optional separate transparent legend
  2) Boxplot panel only (legend removed)
  3) CSV table with disease proximity scores

Defaults target the leo environment layout.
"""
import os
import os
import sys
import argparse
from typing import Dict, Optional, Tuple
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
try:
    import networkx as nx
except ImportError:
    nx = None
DEFAULT_EMBED_PATH = '${SPOKE_DATA_ROOT}/premature_menopause_embedding.csv'
DEFAULT_OUTDIR = '${SPOKE_DATA_ROOT}/biorouter/figs'
DEFAULT_PM_ID = 'Disease:DOID:10787'
DEFAULT_DISEASE_MAP = '${SPOKE_DATA_ROOT}/pm_all_diseases_proximity.csv'
CATEGORIES = ['Other Diseases', 'Cancer', 'Cardiovascular', 'Metabolic', 'Reproductive/Endocrine', 'Bone Diseases', 'OSTEOPOROSIS']
COLORS = {'Other Diseases': '#B0BEC5', 'Cancer': '#E64A19', 'Cardiovascular': '#1E88E5', 'Metabolic': '#8E24AA', 'Reproductive/Endocrine': '#43A047', 'Bone Diseases': '#6D4C41', 'OSTEOPOROSIS': '#D32F2F'}
CANDIDATE_GRAPHS = ['${SPOKE_DATA_ROOT}/data/graph.gpickle', '${SPOKE_DATA_ROOT}/braianpeetoom@Braians-MacBook-Pro/spoke_35M.gpickle', '${SPOKE_DATA_ROOT}/V1_homosapiens_1N_9E.gpickle', '${SPOKE_DATA_ROOT}/extract_spoke/spokeML/extract_spoke/V6_graph_virus.gpickle', '${SPOKE_DATA_ROOT}/ViroLink/VirusLink_project/current_metagraph.gpickle']

def ensure_outdir(path: str):
    os.makedirs(path, exist_ok=True)

def locate_graph(custom_path: Optional[str]=None) -> Optional[str]:
    if custom_path:
        return custom_path if os.path.exists(custom_path) else None
    for p in CANDIDATE_GRAPHS:
        if os.path.exists(p):
            return p
    return None

def load_graph(graph_path: str):
    if nx is None:
        raise RuntimeError('networkx not available; cannot load graph to compute embedding')
    import pickle
    with open(graph_path, 'rb') as f:
        G = pickle.load(f)
    return G

def compute_pm_embedding(G, pm_id: str, alpha: float=0.85) -> pd.DataFrame:
    """Compute personalized PageRank embedding from PM node."""
    if nx is None:
        raise RuntimeError('networkx required to compute embedding')
    if pm_id not in G:
        raise ValueError(f'PM id {pm_id} not found in graph')
    pr = nx.pagerank(G, alpha=alpha, personalization={pm_id: 1})
    df = pd.DataFrame({'node': list(pr.keys()), 'pagerank': list(pr.values())})
    return df

def normalize_doid(node_id: str) -> Optional[str]:
    if node_id is None:
        return None
    s = str(node_id)
    if s.startswith('Disease:DOID:'):
        return s.split('Disease:', 1)[1]
    if s.startswith('DOID:'):
        return s
    return None

def extract_disease_nodes(G) -> Tuple[Dict[str, str], set]:
    """Return mapping node->name and set of disease nodes using graph attributes."""
    name_map = {}
    disease_nodes = set()
    for n, d in G.nodes(data=True):
        nm = d.get('name') or d.get('label') or str(n)
        name_map[n] = nm
        if d.get('node_type') == 'Disease' or (isinstance(n, str) and (n.startswith('Disease:') or n.startswith('DOID:'))):
            disease_nodes.add(n)
    return (name_map, disease_nodes)

def load_doid_name_map(path: Optional[str]) -> Dict[str, str]:
    mapping = {}
    if path and os.path.exists(path):
        try:
            df = pd.read_csv(path)
            id_col = None
            name_col = None
            for c in df.columns:
                cl = c.lower()
                if cl in ('disease_id', 'doid', 'id'):
                    id_col = c
                if cl in ('disease', 'name', 'disease_name'):
                    name_col = c
            if id_col and name_col:
                for _, row in df.iterrows():
                    did = str(row[id_col])
                    nm = str(row[name_col])
                    if did.startswith('DOID:'):
                        mapping[did] = nm
        except Exception as e:
            print(f'WARNING: failed to load disease map from {path}: {e}')
    return mapping

def classify_disease(name: str) -> str:
    low = (name or '').lower()
    if 'dementia' in low:
        return 'DEMENTIA'
    if 'osteoporosis' in low:
        return 'OSTEOPOROSIS'
    if any((k in low for k in ['fracture', 'osteopenia', 'bone', 'skeletal'])):
        return 'Bone Diseases'
    if any((k in low for k in ['cancer', 'neoplasm', 'carcinoma', 'tumor', 'malign'])):
        return 'Cancer'
    if any((k in low for k in ['cardio', 'myocard', 'vascular', 'stroke', 'hypertens', 'arrhythm', 'coronary'])):
        return 'Cardiovascular'
    if any((k in low for k in ['diabet', 'metabolic', 'obes', 'lipid', 'insulin', 'dyslipid'])):
        return 'Metabolic'
    if any((k in low for k in ['ovarian', 'menopause', 'estrogen', 'reproductive', 'endocr', 'pregnan', 'fertil', 'uter', 'mammary'])):
        return 'Reproductive/Endocrine'
    return 'Other Diseases'

def prepare_disease_scores(embedding_df: pd.DataFrame, name_map: Dict[str, str], disease_nodes: set, pm_id: str) -> pd.DataFrame:
    df = embedding_df.copy()
    mask = df['node'].astype(str).str.startswith('DOID:') | df['node'].astype(str).str.startswith('Disease:DOID:')
    df = df[mask].copy()
    df = df[df['node'] != pm_id]
    df['name'] = [name_map.get(normalize_doid(n) or n, str(n)) for n in df['node'].astype(str)]
    df = df.sort_values('pagerank', ascending=False).reset_index(drop=True)
    df['rank'] = np.arange(1, len(df) + 1)
    df['percentile'] = df['pagerank'].rank(pct=True, ascending=True) * 100.0
    df['log_pagerank'] = np.log10(df['pagerank'] + 1e-10)
    df['category'] = [classify_disease(nm) for nm in df['name'].tolist()]
    return df

def find_osteoporosis_row(df: pd.DataFrame) -> Optional[pd.Series]:
    o = df[df['category'] == 'OSTEOPOROSIS']
    if len(o) == 0:
        o = df[df['name'].str.lower().str.contains('osteoporosis', na=False)]
    if len(o) == 0:
        return None
    return o.iloc[0]

def find_dementia_row(df: pd.DataFrame) -> Optional[pd.Series]:
    d = df[df['node'].isin(['DOID:1307', 'Disease:DOID:1307']) | df['name'].str.lower().str.contains('dementia', na=False)]
    if len(d) == 0:
        return None
    return d.iloc[0]

def save_transparent_legend(handles, labels, out_path, ncol=1, fontsize=10, title=None):
    fig = plt.figure(figsize=(6, 0.5 * max(1, len(labels) / ncol)))
    fig.patch.set_alpha(0.0)
    ax = fig.add_subplot(111)
    ax.axis('off')
    ax.legend(handles=handles, labels=labels, loc='center', ncol=ncol, frameon=False, fontsize=fontsize, title=title)
    fig.savefig(out_path, dpi=600, bbox_inches='tight', transparent=True)
    plt.close(fig)

def create_hist_only(disease_scores: pd.DataFrame, osteo_row: Optional[pd.Series], out_png: Optional[str]=None, out_legend: Optional[str]=None, dpi: int=600, normalize_to_one: bool=True):
    plt.rcParams.update({'figure.facecolor': 'white', 'savefig.facecolor': 'white'})
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    trimmed = disease_scores[np.isfinite(disease_scores['log_pagerank'])].copy()
    legend_handles = []
    legend_labels = []
    bins = 40
    for cat in CATEGORIES:
        vals = trimmed[trimmed['category'] == cat]['log_pagerank'].values
        if len(vals) == 0:
            continue
        n, b, patches = ax.hist(vals, bins=bins, alpha=0.6, color=COLORS.get(cat, '#999999'), edgecolor='#FFFFFF', linewidth=0.6, density=True)
        if normalize_to_one and np.max(n) > 0:
            scale = float(np.max(n))
            for p in patches:
                p.set_height(p.get_height() / scale)
        legend_handles.append(plt.Rectangle((0, 0), 1, 1, facecolor=COLORS.get(cat, '#999999'), alpha=0.6))
        legend_labels.append(cat)
    if osteo_row is not None:
        rank_val = int(osteo_row['rank']) if 'rank' in osteo_row else int(getattr(osteo_row, 'rank', -1))
        x = float(osteo_row['log_pagerank'])
        ax.axvline(x, color=COLORS['OSTEOPOROSIS'], linestyle='--', linewidth=2.5)
        legend_handles.append(plt.Line2D([0], [0], color=COLORS['OSTEOPOROSIS'], linestyle='--', linewidth=2.5))
        legend_labels.append('Osteoporosis')
    dementia_row = find_dementia_row(disease_scores)
    if dementia_row is not None:
        try:
            d_rank = int(dementia_row['rank']) if 'rank' in dementia_row else int(getattr(dementia_row, 'rank', -1))
        except Exception:
            d_rank = -1
        dx = float(dementia_row.get('log_pagerank', np.log10(float(dementia_row['pagerank']) + 1e-10)))
        dementia_color = '#FF6F00'
        ax.axvline(dx, color=dementia_color, linestyle='--', linewidth=2.5)
        legend_handles.append(plt.Line2D([0], [0], color=dementia_color, linestyle='--', linewidth=2.5))
        legend_labels.append('Dementia')
    ax.set_xlabel('log₁₀(PageRank + 1e-10)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Density', fontsize=11, fontweight='bold')
    ax.set_title('Disease Network Proximity to Premature Menopause\n(PageRank embedding from PM node)', fontsize=13, fontweight='bold')
    ax.set_xlim(-8, -2)
    ax.grid(True, alpha=0.25)
    if normalize_to_one:
        ax.set_ylim(0, 1.0)
    if out_png:
        fig.savefig(out_png, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    if out_legend:
        save_transparent_legend(legend_handles, legend_labels, out_legend, ncol=1, fontsize=10)

def create_box_only(disease_scores: pd.DataFrame, out_png: Optional[str]=None, dpi: int=600):
    plt.rcParams.update({'figure.facecolor': 'white', 'savefig.facecolor': 'white'})
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    trimmed = disease_scores[np.isfinite(disease_scores['log_pagerank'])].copy()
    box_data, box_labels, box_colors = ([], [], [])
    for cat in CATEGORIES:
        vals = trimmed[trimmed['category'] == cat]['log_pagerank']
        if len(vals) > 0:
            box_data.append(vals)
            box_labels.append(f'{cat}\n(n={len(vals)})')
            box_colors.append(COLORS.get(cat, '#999999'))
    bp = ax.boxplot(box_data, tick_labels=box_labels, patch_artist=True, widths=0.6, showfliers=False)
    for patch, color in zip(bp['boxes'], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
        patch.set_edgecolor('#333333')
    for median in bp['medians']:
        median.set_color('#333333')
        median.set_linewidth(1.2)
    ax.set_ylabel('log₁₀(PageRank + 1e-10)', fontsize=11, fontweight='bold')
    ax.set_title('Proximity by Disease Category', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.25)
    if out_png:
        fig.savefig(out_png, dpi=dpi, bbox_inches='tight')
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--embedding', default=DEFAULT_EMBED_PATH, help='Path to PM-seeded PageRank embedding CSV')
    ap.add_argument('--graph', default=None, help='Path to SPOKE gpickle (optional; required if --compute)')
    ap.add_argument('--outdir', default=DEFAULT_OUTDIR, help='Directory to save outputs')
    ap.add_argument('--pm-id', default=DEFAULT_PM_ID, help='Premature Menopause node id in graph')
    ap.add_argument('--dpi', type=int, default=600, help='DPI for saved figure')
    ap.add_argument('--compute', action='store_true', help='Compute embedding from graph instead of loading CSV')
    ap.add_argument('--alpha', type=float, default=0.85, help='PageRank alpha (if computing)')
    ap.add_argument('--skip-graph', dest='skip_graph', action='store_true', help='Skip loading graph for disease names; use DOID->name map instead')
    ap.add_argument('--disease-map', default=DEFAULT_DISEASE_MAP, help='CSV with DOID->disease name (e.g., pm_all_diseases_proximity.csv)')
    ap.add_argument('--out-hist', dest='out_hist', default=None, help='Path to save histogram-only panel (PNG)')
    ap.add_argument('--out-box', dest='out_box', default=None, help='Path to save boxplot-only panel (PNG)')
    ap.add_argument('--out-legend', dest='out_legend', default=None, help='Path to save legend (transparent PNG)')
    args = ap.parse_args()
    ensure_outdir(args.outdir)
    if args.compute:
        gp = locate_graph(args.graph)
        if gp is None:
            print('ERROR: No graph path found; provide --graph')
            sys.exit(1)
        print(f'Loading graph: {gp}')
        G = load_graph(gp)
        print(f'Graph nodes: {G.number_of_nodes():,}; edges: {G.number_of_edges():,}')
        print(f'Computing personalized PageRank from {args.pm_id} (alpha={args.alpha})…')
        embedding_df = compute_pm_embedding(G, args.pm_id, alpha=args.alpha)
        name_map, disease_nodes = extract_disease_nodes(G)
    else:
        if not os.path.exists(args.embedding):
            print(f'ERROR: Embedding not found: {args.embedding}. Use --compute with --graph.')
            sys.exit(1)
        print(f'Loading embedding: {args.embedding}')
        embedding_df = pd.read_csv(args.embedding)
        if set(embedding_df.columns) != {'node', 'pagerank'}:
            for cnode in ['node', 'id', 'name']:
                for cpr in ['pagerank', 'pr', 'score', 'ppr']:
                    if cnode in embedding_df.columns and cpr in embedding_df.columns:
                        embedding_df = embedding_df[[cnode, cpr]].rename(columns={cnode: 'node', cpr: 'pagerank'})
                        break
                if set(embedding_df.columns) == {'node', 'pagerank'}:
                    break
        name_map = {}
        disease_nodes = set()
        doid_map = load_doid_name_map(args.disease_map)
        if args.skip_graph:
            print('Skipping graph load; using DOID->name map (if available)')
        else:
            gp = locate_graph(args.graph)
            if gp is not None:
                print(f'Loading graph for disease mapping: {gp}')
                G = load_graph(gp)
                name_map_graph, disease_nodes_graph = extract_disease_nodes(G)
                for n, nm in name_map_graph.items():
                    did = normalize_doid(n)
                    if did:
                        name_map[did] = nm
                disease_nodes = set(disease_nodes_graph)
            else:
                print('Graph not found. Falling back to DOID->name map only.')
        name_map.update(doid_map)
        if not disease_nodes:
            disease_nodes = set([n for n in embedding_df['node'].astype(str) if n.startswith('DOID:') or n.startswith('Disease:DOID:')])
    disease_scores = prepare_disease_scores(embedding_df, name_map, disease_nodes, pm_id=args.pm_id)
    out_csv = os.path.join(args.outdir, 'pm_disease_proximity_scores.csv')
    disease_scores.to_csv(out_csv, index=False)
    print(f'Saved: {out_csv} (rows: {len(disease_scores):,})')
    out_hist = args.out_hist or os.path.join(args.outdir, 'pm_disease_proximity_hist.png')
    out_box = args.out_box or os.path.join(args.outdir, 'pm_disease_proximity_box.png')
    out_legend = args.out_legend or os.path.join(args.outdir, 'pm_disease_proximity_legend.png')
    create_hist_only(disease_scores, find_osteoporosis_row(disease_scores), out_png=out_hist, out_legend=out_legend, dpi=args.dpi, normalize_to_one=True)
    print(f'Saved: {out_hist}')
    print(f'Saved: {out_legend} (transparent)')
    create_box_only(disease_scores, out_png=out_box, dpi=args.dpi)
    print(f'Saved: {out_box}')
    top15 = disease_scores.head(15)[['rank', 'name', 'category', 'pagerank', 'percentile']]
    print('\nTop 15 diseases by proximity to PM:')
    for _, r in top15.iterrows():
        mark = ' → OSTEOPOROSIS' if r['category'] == 'OSTEOPOROSIS' else ''
        print(f"{int(r['rank']):2d}. {r['name'][:40]:<40} {r['category']:<22} {r['pagerank']:.2e} (Top {r['percentile']:.1f}%)" + mark)
if __name__ == '__main__':
    main()
