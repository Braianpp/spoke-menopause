# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
import pandas as pd
import numpy as np
import networkx as nx
from scipy.stats import fisher_exact
from collections import Counter
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')
try:
    G
except NameError:
    try:
        G = g
        print('Using graph object: g → assigned to G')
    except NameError:
        raise NameError('Could not find graph object. Please load SPOKE as `G` or `g` first.')
TOP_PERCENT = 2
MAX_HOPS = 2
SAVE_PATH = Path('${SPOKE_DATA_ROOT}/')
SAVE_PATH.mkdir(parents=True, exist_ok=True)
suffix = f'_top{TOP_PERCENT}pct_{MAX_HOPS}hops'
print('=' * 80)
print('Configuration')
print('=' * 80)
print(f'Top percent: {TOP_PERCENT}%')
print(f'Max hops: {MAX_HOPS}')
print(f'Save path: {SAVE_PATH}')
print(f'Output suffix: {suffix}')
adverse_pregnancy_outcomes = {'DOID:11714': 'gestational diabetes', 'DOID:10591': 'pre-eclampsia', 'DOID:10590': 'mild pre-eclampsia', 'DOID:13129': 'severe pre-eclampsia', 'DOID:13593': 'eclampsia', 'DOID:13133': 'HELLP syndrome', 'DOID:9667': 'placental abruption'}
menopause_symptoms = {'D019584': 'Hot flashes', 'D006945': 'Hyperhidrosis / sweating', 'D007319': 'Insomnia', 'D012893': 'Sleep-wake disorders', 'D003866': 'Depression', 'D001007': 'Anxiety', 'D008569': 'Memory disorders', 'D060825': 'Cognitive impairment', 'D005221': 'Fatigue', 'D006973': 'Hypertension', 'D009765': 'Obesity', 'D007333': 'Insulin resistance', 'D015430': 'Weight gain', 'D010024': 'Osteoporosis', 'D014549': 'Urinary incontinence', 'D004414': 'Dyspareunia'}
cognitive_targets = {'DOID:10652': "Alzheimer's disease", 'DOID:1307': 'dementia', 'DOID:8725': 'vascular dementia', 'DOID:0080831': 'subjective cognitive decline'}
condition_groups = {'APO': list(adverse_pregnancy_outcomes.keys()), 'Menopause': list(menopause_symptoms.keys()), 'Cognitive': list(cognitive_targets.keys())}
all_labels = {**adverse_pregnancy_outcomes, **menopause_symptoms, **cognitive_targets}
print('\n' + '=' * 80)
print('Concept mapping to SPOKE')
print('=' * 80)
mapped_concepts = []
missing_concepts = []
for group_name, node_ids in condition_groups.items():
    print(f'\n{group_name}')
    print('-' * 80)
    for node_id in node_ids:
        expected_label = all_labels.get(node_id, 'NA')
        if node_id in G.nodes:
            node_data = G.nodes[node_id]
            spoke_name = node_data.get('name', 'NA')
            node_type = node_data.get('node_type', 'NA')
            degree = G.degree(node_id)
            mapped_concepts.append({'group': group_name, 'node_id': node_id, 'expected_label': expected_label, 'spoke_name': spoke_name, 'node_type': node_type, 'degree': degree, 'mapped': True})
            print(f'  ✓ {node_id:15s} | {expected_label:35s} | SPOKE: {spoke_name} [{node_type}] degree={degree}')
        else:
            missing_concepts.append({'group': group_name, 'node_id': node_id, 'expected_label': expected_label, 'spoke_name': None, 'node_type': None, 'degree': None, 'mapped': False})
            print(f'  ✗ {node_id:15s} | {expected_label:35s} | NOT FOUND')
df_concept_mapping = pd.DataFrame(mapped_concepts + missing_concepts)
df_concept_mapping.to_csv(SAVE_PATH / f'concept_mapping_v4_compute{suffix}.csv', index=False)
print('\nMapped concepts:', len(mapped_concepts))
print('Missing concepts:', len(missing_concepts))
print(f'Saved: concept_mapping_v4_compute{suffix}.csv')

def get_connected_genes(G, seed_nodes, max_hops=2):
    """
    BFS from seed nodes and return all Gene nodes within max_hops.
    Returns dictionary: gene_id -> minimum distance from any seed.
    """
    valid_seeds = [n for n in seed_nodes if n in G.nodes]
    connected_genes = {}
    for seed in valid_seeds:
        visited = {seed: 0}
        queue = [(seed, 0)]
        while queue:
            current, depth = queue.pop(0)
            if depth >= max_hops:
                continue
            for neighbor in G.neighbors(current):
                if neighbor in visited:
                    continue
                visited[neighbor] = depth + 1
                queue.append((neighbor, depth + 1))
                if G.nodes[neighbor].get('node_type', '') == 'Gene':
                    if neighbor not in connected_genes:
                        connected_genes[neighbor] = depth + 1
                    else:
                        connected_genes[neighbor] = min(connected_genes[neighbor], depth + 1)
    return connected_genes

def generate_ppr_ranked_genes(G, seed_nodes, connected_genes, alpha=0.85, max_iter=500):
    """
    Run Personalized PageRank over the full graph, then restrict output to connected genes.
    """
    valid_seeds = [n for n in seed_nodes if n in G.nodes]
    if not valid_seeds:
        return pd.DataFrame(columns=['gene_id', 'gene_name', 'pagerank', 'distance', 'rank', 'percentile'])
    personalization = {node: 1.0 for node in valid_seeds}
    pagerank = nx.pagerank(G, alpha=alpha, personalization=personalization, max_iter=max_iter, tol=1e-10)
    rows = []
    for gene_id, distance in connected_genes.items():
        if gene_id not in pagerank:
            continue
        rows.append({'gene_id': gene_id, 'gene_name': G.nodes[gene_id].get('name', gene_id), 'pagerank': pagerank[gene_id], 'distance': distance})
    ranked_df = pd.DataFrame(rows)
    if ranked_df.empty:
        return ranked_df
    ranked_df = ranked_df.sort_values('pagerank', ascending=False).reset_index(drop=True)
    ranked_df['rank'] = np.arange(1, len(ranked_df) + 1)
    ranked_df['percentile'] = (1 - ranked_df['rank'] / len(ranked_df)) * 100
    return ranked_df

def extract_top_genes(ranked_df, top_percent):
    """
    Extract the top X% of a ranked gene dataframe.
    """
    if ranked_df is None or ranked_df.empty:
        return (set(), pd.DataFrame())
    n_top = int(np.ceil(len(ranked_df) * (top_percent / 100.0)))
    top_df = ranked_df.head(n_top).copy()
    gene_set = set(top_df['gene_id'].values)
    return (gene_set, top_df)

def jaccard(s1, s2):
    union = s1 | s2
    return len(s1 & s2) / len(union) if union else 0.0

def gene_set_to_df(gene_set):
    rows = []
    for gene_id in sorted(gene_set):
        gene_name = G.nodes[gene_id].get('name', gene_id) if gene_id in G.nodes else gene_id
        rows.append({'gene_id': gene_id, 'gene_name': gene_name})
    return pd.DataFrame(rows)
print('\n' + '=' * 80)
print('Computing 2-hop connected genes and PPR rankings')
print('=' * 80)
connected_genes_dict = {}
ranked_dfs = {}
gene_sets = {}
top_gene_dfs = {}
for group_name, node_ids in condition_groups.items():
    print(f'\n{group_name}')
    print('-' * 80)
    connected = get_connected_genes(G, node_ids, max_hops=MAX_HOPS)
    connected_genes_dict[group_name] = connected
    dist_counts = Counter(connected.values())
    print(f'Connected genes within {MAX_HOPS} hops: {len(connected):,}')
    for dist in sorted(dist_counts):
        print(f'  {dist} hop(s): {dist_counts[dist]:,}')
    ranked_df = generate_ppr_ranked_genes(G, node_ids, connected)
    ranked_dfs[group_name] = ranked_df
    if ranked_df.empty:
        print('No ranked genes found.')
        gene_sets[group_name] = set()
        top_gene_dfs[group_name] = pd.DataFrame()
        continue
    print(f'PPR-ranked genes: {len(ranked_df):,}')
    print(f"PageRank range: {ranked_df['pagerank'].min():.2e} → {ranked_df['pagerank'].max():.2e}")
    genes, top_df = extract_top_genes(ranked_df, TOP_PERCENT)
    gene_sets[group_name] = genes
    top_gene_dfs[group_name] = top_df
    print(f'Selected top {TOP_PERCENT}%: {len(genes):,} genes')
    print('Top 5 genes:')
    for _, row in top_df.head(5).iterrows():
        print(f"  • {row['gene_name']:<20s} rank={row['rank']:<5d} distance={row['distance']} PR={row['pagerank']:.2e}")
    ranked_df.to_csv(SAVE_PATH / f'{group_name.lower()}_ranked_genes_v4_compute_top100pct_{MAX_HOPS}hops.csv', index=False)
    top_df.to_csv(SAVE_PATH / f'{group_name.lower()}_ranked_genes_v4_compute{suffix}.csv', index=False)
apo_genes = gene_sets['APO']
meno_genes = gene_sets['Menopause']
cog_genes = gene_sets['Cognitive']
apo_meno = apo_genes & meno_genes
apo_cog = apo_genes & cog_genes
meno_cog = meno_genes & cog_genes
all_three = apo_genes & meno_genes & cog_genes
jaccard_matrix = pd.DataFrame([[1.0, jaccard(apo_genes, meno_genes), jaccard(apo_genes, cog_genes)], [jaccard(meno_genes, apo_genes), 1.0, jaccard(meno_genes, cog_genes)], [jaccard(cog_genes, apo_genes), jaccard(cog_genes, meno_genes), 1.0]], index=['APO', 'Menopause', 'Cognitive'], columns=['APO', 'Menopause', 'Cognitive'])
print('\n' + '=' * 80)
print('Overlap analysis')
print('=' * 80)
print(f'APO genes:        {len(apo_genes):,}')
print(f'Menopause genes:  {len(meno_genes):,}')
print(f'Cognitive genes:  {len(cog_genes):,}')
print(f'\nAPO ∩ Menopause:     {len(apo_meno):,} | Jaccard={jaccard(apo_genes, meno_genes):.4f}')
print(f'APO ∩ Cognitive:     {len(apo_cog):,} | Jaccard={jaccard(apo_genes, cog_genes):.4f}')
print(f'Menopause ∩ Cognitive: {len(meno_cog):,} | Jaccard={jaccard(meno_genes, cog_genes):.4f}')
print(f'All three:           {len(all_three):,}')
if all_three:
    print('\nTop shared genes across all three:')
    shared_names = gene_set_to_df(all_three).sort_values('gene_name')
    for _, row in shared_names.head(30).iterrows():
        print(f"  • {row['gene_name']} ({row['gene_id']})")
print('\n' + '=' * 80)
print('Fisher exact tests')
print('=' * 80)
all_graph_genes = set((n for n in G.nodes if G.nodes[n].get('node_type') == 'Gene'))
print(f'Background genes in SPOKE: {len(all_graph_genes):,}')
fisher_rows = []

def fisher_enrichment(set_a, set_b, background, label):
    a = len(set_a & set_b)
    b = len(set_a - set_b)
    c = len(set_b - set_a)
    d = len(background - set_a - set_b)
    table = np.array([[a, b], [c, d]])
    odds_ratio, p_value = fisher_exact(table, alternative='greater')
    return {'test': label, 'a_overlap': a, 'b_a_only': b, 'c_b_only': c, 'd_neither': d, 'odds_ratio': odds_ratio, 'p_value': p_value}
fisher_rows.append(fisher_enrichment(apo_genes, meno_cog, all_graph_genes, 'APO enriched in Menopause ∩ Cognitive'))
fisher_rows.append(fisher_enrichment(apo_genes, meno_genes, all_graph_genes, 'APO enriched in Menopause'))
fisher_rows.append(fisher_enrichment(apo_genes, cog_genes, all_graph_genes, 'APO enriched in Cognitive'))
fisher_results = pd.DataFrame(fisher_rows)
for _, row in fisher_results.iterrows():
    sig = '***' if row['p_value'] < 0.001 else '**' if row['p_value'] < 0.01 else '*' if row['p_value'] < 0.05 else 'ns'
    print(f"{row['test']}: OR={row['odds_ratio']:.2f}, p={row['p_value']:.2e} {sig}")
print('\n' + '=' * 80)
print('Sensitivity sweep across PPR thresholds')
print('=' * 80)
sweep_percents = [1, 2, 5, 10, 15, 20, 30, 50, 75, 100]
sweep_rows = []
for pct in sweep_percents:
    sweep_sets = {}
    for group_name, ranked_df in ranked_dfs.items():
        genes_pct, _ = extract_top_genes(ranked_df, pct)
        sweep_sets[group_name] = genes_pct
    am = sweep_sets['APO'] & sweep_sets['Menopause']
    ac = sweep_sets['APO'] & sweep_sets['Cognitive']
    mc = sweep_sets['Menopause'] & sweep_sets['Cognitive']
    a3 = sweep_sets['APO'] & sweep_sets['Menopause'] & sweep_sets['Cognitive']
    sweep_rows.append({'top_percent': pct, 'apo_genes': len(sweep_sets['APO']), 'menopause_genes': len(sweep_sets['Menopause']), 'cognitive_genes': len(sweep_sets['Cognitive']), 'apo_meno': len(am), 'apo_cog': len(ac), 'meno_cog': len(mc), 'all_three': len(a3), 'jaccard_apo_meno': jaccard(sweep_sets['APO'], sweep_sets['Menopause']), 'jaccard_apo_cog': jaccard(sweep_sets['APO'], sweep_sets['Cognitive']), 'jaccard_meno_cog': jaccard(sweep_sets['Menopause'], sweep_sets['Cognitive'])})
df_sweep = pd.DataFrame(sweep_rows)
print(df_sweep.to_string(index=False))
print('\n' + '=' * 80)
print('Saving outputs')
print('=' * 80)
gene_set_to_df(apo_genes).to_csv(SAVE_PATH / f'apo_genes_v4_compute{suffix}.csv', index=False)
gene_set_to_df(meno_genes).to_csv(SAVE_PATH / f'menopause_genes_v4_compute{suffix}.csv', index=False)
gene_set_to_df(cog_genes).to_csv(SAVE_PATH / f'cognitive_genes_v4_compute{suffix}.csv', index=False)
gene_set_to_df(apo_meno).to_csv(SAVE_PATH / f'shared_apo_meno_v4_compute{suffix}.csv', index=False)
gene_set_to_df(apo_cog).to_csv(SAVE_PATH / f'shared_apo_cog_v4_compute{suffix}.csv', index=False)
gene_set_to_df(meno_cog).to_csv(SAVE_PATH / f'shared_meno_cog_v4_compute{suffix}.csv', index=False)
gene_set_to_df(all_three).to_csv(SAVE_PATH / f'shared_all_three_v4_compute{suffix}.csv', index=False)
fisher_results.to_csv(SAVE_PATH / f'fisher_exact_results_v4_compute{suffix}.csv', index=False)
df_sweep.to_csv(SAVE_PATH / f'sensitivity_sweep_v4_compute_{MAX_HOPS}hops.csv', index=False)
jaccard_matrix.to_csv(SAVE_PATH / f'jaccard_matrix_v4_compute{suffix}.csv')
summary = pd.DataFrame({'Group': ['APO', 'Menopause', 'Cognitive', 'APO∩Menopause', 'APO∩Cognitive', 'Menopause∩Cognitive', 'All three'], 'Total_connected': [len(connected_genes_dict['APO']), len(connected_genes_dict['Menopause']), len(connected_genes_dict['Cognitive']), np.nan, np.nan, np.nan, np.nan], 'Selected_genes': [len(apo_genes), len(meno_genes), len(cog_genes), len(apo_meno), len(apo_cog), len(meno_cog), len(all_three)], 'Top_percent': TOP_PERCENT, 'Max_hops': MAX_HOPS})
summary.to_csv(SAVE_PATH / f'summary_v4_compute{suffix}.csv', index=False)
print(f'Saved all outputs to: {SAVE_PATH}')
print('\nReady for plotting cell.')
