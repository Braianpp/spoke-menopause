# Credential-free source extraction. See docs/SOURCE_MAP.csv.
"""
Pharmacogene × MHT Compound Interaction Map from SPOKE
======================================================
Visualizes: Which pharmacogenes from PharmGKB encode proteins
that physically bind or interact with specific MHT compounds?

Requires: access to SPOKE via Cypher (open_spoke) or pre-loaded G
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm
import seaborn as sns
PHARMGKB_GENES = ['ESR1', 'CYP2B6', 'CYP2C19', 'CYP2C9', 'SULT1A1', 'UGT1A1', 'VDR', 'TP53', 'NFE2L2', 'NAT2', 'TCL1A', 'ABCC1', 'SLCO1B1', 'ADAMTS10']
MHT_COMPOUNDS = {'inchikey:VOXZDWNPVJITMN-ZBRFXRBCSA-N': {'name': 'Estradiol', 'class': 'Estrogen'}, 'inchikey:FHXBMXJMKMWVRG-SLHNCBLASA-N': {'name': 'Estradiol Acetate', 'class': 'Estrogen'}, 'inchikey:BFPYWIDHMRZLRN-SLHNCBLASA-N': {'name': 'Ethinyl Estradiol', 'class': 'Estrogen'}, 'inchikey:DNXHEGUUPJUMQT-CBZIJGRNSA-N': {'name': 'Estrone', 'class': 'Estrogen (CE)'}, 'inchikey:WKRLQDKEXYKHJB-HFTRVMKXSA-N': {'name': 'Equilin', 'class': 'Estrogen (CE)'}, 'inchikey:PDRGHUMCVRDZLQ-WMZOPIPTSA-N': {'name': 'Equilenin', 'class': 'Estrogen (CE)'}, 'inchikey:UCJGJABZCDBEDK-UHFFFAOYSA-N': {'name': 'Bazedoxifene', 'class': 'SERM'}}
MHT_IDS = list(MHT_COMPOUNDS.keys())
GENE_CATEGORIES = {'ESR1': 'Receptor', 'VDR': 'Receptor', 'CYP2B6': 'Phase I metabolism', 'CYP2C19': 'Phase I metabolism', 'CYP2C9': 'Phase I metabolism', 'SULT1A1': 'Phase II metabolism', 'UGT1A1': 'Phase II metabolism', 'NAT2': 'Phase II metabolism', 'ABCC1': 'Transporter', 'SLCO1B1': 'Transporter', 'TP53': 'Tumor suppressor', 'NFE2L2': 'Stress response', 'TCL1A': 'Signaling', 'ADAMTS10': 'Other'}
GENE_DESCRIPTIONS = {'ESR1': 'Estrogen receptor alpha — primary target of estrogen MHT', 'VDR': 'Vitamin D receptor — nuclear receptor, bone metabolism', 'CYP2B6': 'Cytochrome P450 2B6 — oxidative metabolism of estradiol', 'CYP2C19': 'Cytochrome P450 2C19 — metabolism of norethindrone, progesterone', 'CYP2C9': 'Cytochrome P450 2C9 — metabolism of progesterone, estradiol', 'SULT1A1': 'Sulfotransferase 1A1 — sulfation of estradiol, ethinyl estradiol', 'UGT1A1': 'UDP-glucuronosyltransferase 1A1 — glucuronidation of ethinyl estradiol, estrone', 'NAT2': 'N-acetyltransferase 2 — phase II acetylation enzyme', 'ABCC1': 'ATP-binding cassette transporter C1 — drug efflux', 'SLCO1B1': 'Solute carrier organic anion transporter 1B1 — hepatic uptake', 'TP53': 'Tumor protein p53 — tumor suppressor, cell cycle regulation', 'NFE2L2': 'Nuclear factor erythroid 2-related factor 2 — oxidative stress response', 'TCL1A': 'T-cell leukemia/lymphoma protein 1A — estrogen-responsive signaling', 'ADAMTS10': 'Metalloproteinase with thrombospondin motifs 10'}
CAT_ORDER = ['Receptor', 'Phase I metabolism', 'Phase II metabolism', 'Transporter', 'Tumor suppressor', 'Stress response', 'Signaling', 'Other']

def query_spoke_pharmacogene_mht(open_spoke, SPOKE_User, SPOKE_Password, SPOKE_DEV_URI):
    mht_ids_str = '[' + ','.join([f"'{x}'" for x in MHT_IDS]) + ']'
    gene_str = '[' + ','.join([f"'{g}'" for g in PHARMGKB_GENES]) + ']'
    q_binds = f"\n    MATCH (c:Compound)-[r:BINDS_CbP]->(p:Protein)<-[:ENCODES_GeP]-(g:Gene)\n    WHERE c.identifier IN {mht_ids_str}\n      AND g.name IN {gene_str}\n    RETURN DISTINCT c.identifier AS compound_id, c.name AS compound,\n           p.name AS protein, g.name AS gene, 'BINDS' AS edge_type\n    "
    q_interacts = f"\n    MATCH (p:Protein)-[r:INTERACTS_PiC]->(c:Compound), (g:Gene)-[:ENCODES_GeP]->(p)\n    WHERE c.identifier IN {mht_ids_str}\n      AND g.name IN {gene_str}\n    RETURN DISTINCT c.identifier AS compound_id, c.name AS compound,\n           p.name AS protein, g.name AS gene, 'INTERACTS' AS edge_type\n    "
    q_regulates = f"\n    MATCH (c:Compound)-[r]->(g:Gene)\n    WHERE c.identifier IN {mht_ids_str}\n      AND g.name IN {gene_str}\n      AND type(r) IN ['UPREGULATES_CuG', 'DOWNREGULATES_CdG']\n    RETURN DISTINCT c.identifier AS compound_id, c.name AS compound,\n           g.name AS gene,\n           CASE WHEN type(r) = 'UPREGULATES_CuG' THEN 'UPREGULATES'\n                ELSE 'DOWNREGULATES' END AS edge_type\n    "
    dfs = []
    with open_spoke(SPOKE_User, SPOKE_Password, SPOKE_DEV_URI) as db:
        session = db.session()
        for q in [q_binds, q_interacts, q_regulates]:
            results = session.run(q)
            dfs.append(pd.DataFrame(results.data()))
    all_edges = pd.concat(dfs, ignore_index=True)
    all_edges['compound_short'] = all_edges['compound_id'].map(lambda x: MHT_COMPOUNDS.get(x, {}).get('name', x))
    all_edges['compound_class'] = all_edges['compound_id'].map(lambda x: MHT_COMPOUNDS.get(x, {}).get('class', 'Unknown'))
    print(f'✓ Total edges: {len(all_edges)}')
    print(f'  BINDS: {len(dfs[0])}, INTERACTS: {len(dfs[1])}, REGULATES: {len(dfs[2])}')
    print(f"  Unique genes: {all_edges['gene'].nunique()}")
    print(f"  Unique compounds: {all_edges['compound_short'].nunique()}")
    return all_edges

def build_interaction_matrix(all_edges):
    grouped = all_edges.groupby(['gene', 'compound_short'])['edge_type'].apply(set).reset_index()

    def encode(etypes):
        has_bind = 'BINDS' in etypes
        has_inter = 'INTERACTS' in etypes
        has_up = 'UPREGULATES' in etypes
        has_down = 'DOWNREGULATES' in etypes
        if has_bind and has_inter and (has_up or has_down):
            return 7
        if has_bind and has_inter:
            return 3
        if has_bind and (has_up or has_down):
            return 6
        if has_inter and (has_up or has_down):
            return 6
        if has_bind:
            return 1
        if has_inter:
            return 2
        if has_up:
            return 4
        if has_down:
            return 5
        return 0
    grouped['code'] = grouped['edge_type'].apply(encode)
    class_order = ['Estrogen', 'Estrogen (CE)', 'Progestin', 'SERM']
    compound_order = []
    for cls in class_order:
        compound_order.extend(sorted([v['name'] for v in MHT_COMPOUNDS.values() if v['class'] == cls]))
    compounds_with_edges = set(all_edges['compound_short'].unique())
    compound_order = [c for c in compound_order if c in compounds_with_edges]
    genes_with_edges = set(all_edges['gene'].unique())
    gene_order = []
    for cat in CAT_ORDER:
        gene_order.extend(sorted([g for g, c in GENE_CATEGORIES.items() if c == cat and g in genes_with_edges]))
    matrix = pd.DataFrame(0, index=gene_order, columns=compound_order)
    for _, row in grouped.iterrows():
        if row['gene'] in matrix.index and row['compound_short'] in matrix.columns:
            matrix.loc[row['gene'], row['compound_short']] = row['code']
    return matrix

def plot_heatmap(matrix, save_path=''):
    fig, ax = plt.subplots(figsize=(14, 6.5))
    cmap_colors = ['#F8F9FA', '#3B82F6', '#10B981', '#8B5CF6', '#F97316', '#EF4444', '#1E3A5F', '#111827']
    cmap = ListedColormap(cmap_colors)
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]
    norm = BoundaryNorm(bounds, cmap.N)
    im = ax.imshow(matrix.values, cmap=cmap, norm=norm, aspect='auto')
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, rotation=45, ha='right', fontsize=16)
    ax.set_yticks(range(len(matrix.index)))
    ax.set_yticklabels(matrix.index, fontsize=16)
    compound_class_colors = {'Estrogen': '#E74C8B', 'Estrogen (CE)': '#F4A6C9', 'Progestin': '#6366F1', 'SERM': '#F59E0B'}
    for j, comp in enumerate(matrix.columns):
        cls = next((v['class'] for v in MHT_COMPOUNDS.values() if v['name'] == comp), '')
        color = compound_class_colors.get(cls, '#D1D5DB')
        rect = plt.Rectangle((j - 0.4, -1.4), 0.8, 0.6, facecolor=color, edgecolor='none', clip_on=False)
        ax.add_patch(rect)
    for i in range(len(matrix.index) + 1):
        ax.axhline(i - 0.5, color='white', linewidth=1.5)
    for j in range(len(matrix.columns) + 1):
        ax.axvline(j - 0.5, color='white', linewidth=1.5)
    legend_elements = [mpatches.Patch(facecolor='#3B82F6', label='Binds protein'), mpatches.Patch(facecolor='#10B981', label='Interacts with protein'), mpatches.Patch(facecolor='#8B5CF6', label='Binds + Interacts'), mpatches.Patch(facecolor='#F97316', label='Upregulates gene'), mpatches.Patch(facecolor='#EF4444', label='Downregulates gene'), mpatches.Patch(facecolor='#1E3A5F', label='Binding + Regulation')]
    leg1 = ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1.0), fontsize=9, title='Edge Type', title_fontsize=10, frameon=True)
    ax.add_artist(leg1)
    comp_legend = [mpatches.Patch(facecolor=c, label=k) for k, c in compound_class_colors.items()]
    ax.legend(handles=comp_legend, loc='upper left', bbox_to_anchor=(1.02, 0.3), fontsize=9, title='Drug Class', title_fontsize=10, frameon=True)
    ax.set_title('', fontsize=13, fontweight='bold', pad=25)
    plt.tight_layout()
    if save_path:
        fig.savefig(f'{save_path}pharmacogene_mht_heatmap.pdf', dpi=300, bbox_inches='tight')
        fig.savefig(f'{save_path}pharmacogene_mht_heatmap.png', dpi=300, bbox_inches='tight')
        print(f'✓ Heatmap saved to {save_path}')
    return fig

def plot_gene_table(all_edges, save_path=''):
    """
    Clean table figure: gene, function, description, MHT connections, edge types.
    Includes genes with zero SPOKE edges (from PharmGKB but absent in graph).
    """
    genes_with_edges = set(all_edges['gene'].unique())
    rows = []
    for cat in CAT_ORDER:
        genes_in_cat = sorted([g for g, c in GENE_CATEGORIES.items() if c == cat and g in genes_with_edges])
        for gene in genes_in_cat:
            n_compounds = all_edges[all_edges['gene'] == gene]['compound_short'].nunique()
            edge_types = ', '.join(sorted(all_edges[all_edges['gene'] == gene]['edge_type'].unique()))
            rows.append({'Gene': gene, 'Function': cat, 'Description': GENE_DESCRIPTIONS.get(gene, ''), 'MHT connections': n_compounds, 'Edge types': edge_types})
    for gene in PHARMGKB_GENES:
        if gene not in genes_with_edges:
            rows.append({'Gene': gene, 'Function': GENE_CATEGORIES.get(gene, 'Unknown'), 'Description': GENE_DESCRIPTIONS.get(gene, ''), 'MHT connections': 0, 'Edge types': '— not connected in SPOKE'})
    table_df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(16, 0.45 * len(rows) + 1.8))
    ax.axis('off')
    cat_bg = {'Receptor': '#FDEEF4', 'Phase I metabolism': '#EFF6FF', 'Phase II metabolism': '#ECFDF5', 'Transporter': '#FFFBEB', 'Tumor suppressor': '#F5F3FF', 'Stress response': '#FEF2F2', 'Signaling': '#F3F4F6', 'Other': '#F9FAFB'}
    cat_fg = {'Receptor': '#BE185D', 'Phase I metabolism': '#1D4ED8', 'Phase II metabolism': '#059669', 'Transporter': '#D97706', 'Tumor suppressor': '#7C3AED', 'Stress response': '#DC2626', 'Signaling': '#4B5563', 'Other': '#9CA3AF'}
    col_labels = ['Gene', 'Function', 'Description', 'MHT\nConnections', 'Edge Types']
    cell_text = []
    cell_colours = []
    for _, row in table_df.iterrows():
        bg = cat_bg.get(row['Function'], '#FFFFFF')
        cell_text.append([row['Gene'], row['Function'], row['Description'], str(row['MHT connections']), row['Edge types']])
        cell_colours.append([bg] * 5)
    table = ax.table(cellText=cell_text, colLabels=col_labels, cellColours=cell_colours, colColours=['#1E2761'] * 5, cellLoc='left', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.6)
    for j in range(5):
        cell = table[0, j]
        cell.set_text_props(fontweight='bold', color='white', fontsize=10)
        cell.set_edgecolor('#1E2761')
    for i in range(len(cell_text)):
        cat = table_df.iloc[i]['Function']
        for j in range(5):
            cell = table[i + 1, j]
            cell.set_edgecolor('#E5E7EB')
            if j == 0:
                cell.set_text_props(fontweight='bold', fontfamily='monospace', fontsize=10)
            elif j == 1:
                cell.set_text_props(color=cat_fg.get(cat, '#333'), fontweight='bold')
            elif j == 3:
                cell.set_text_props(ha='center')
                cell._loc = 'center'
            if table_df.iloc[i]['MHT connections'] == 0:
                cell.set_text_props(color='#9CA3AF')
    col_widths = [0.07, 0.12, 0.45, 0.08, 0.2]
    for j, w in enumerate(col_widths):
        for i in range(len(cell_text) + 1):
            table[i, j].set_width(w)
    ax.set_title('PharmGKB Pharmacogenes — Functional Classification & SPOKE Connectivity', fontsize=13, fontweight='bold', pad=15)
    plt.tight_layout()
    if save_path:
        fig.savefig(f'{save_path}pharmacogene_table.pdf', dpi=300, bbox_inches='tight')
        fig.savefig(f'{save_path}pharmacogene_table.png', dpi=300, bbox_inches='tight')
        print(f'✓ Gene table saved to {save_path}')
    return (fig, table_df)

def plot_compound_connectivity(all_edges, save_path=''):
    fig, ax = plt.subplots(figsize=(12, 6))
    compound_class_colors = {'Estrogen': '#E74C8B', 'Estrogen (CE)': '#F4A6C9', 'Progestin': '#6366F1', 'SERM': '#F59E0B'}
    counts = all_edges.groupby(['compound_short', 'compound_class', 'edge_type'])['gene'].nunique().reset_index()
    counts.columns = ['compound', 'compound_class', 'edge_type', 'n_genes']
    pivot = counts.pivot_table(index=['compound', 'compound_class'], columns='edge_type', values='n_genes', fill_value=0)
    pivot['total'] = pivot.sum(axis=1)
    pivot = pivot.sort_values('total', ascending=True)
    edge_colors = {'BINDS': '#3B82F6', 'INTERACTS': '#10B981', 'UPREGULATES': '#F97316', 'DOWNREGULATES': '#EF4444'}
    y_pos = range(len(pivot))
    left = np.zeros(len(pivot))
    for etype in ['BINDS', 'INTERACTS', 'UPREGULATES', 'DOWNREGULATES']:
        if etype in pivot.columns:
            vals = pivot[etype].values
            ax.barh(y_pos, vals, left=left, color=edge_colors[etype], label=etype.capitalize(), edgecolor='white', linewidth=0.5, height=0.7)
            left += vals
    for i, ((comp, cls), _) in enumerate(pivot.iterrows()):
        color = compound_class_colors.get(cls, '#D1D5DB')
        ax.plot(-0.3, i, 's', color=color, markersize=10, clip_on=False)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([idx[0] for idx in pivot.index], fontsize=10, fontweight='bold')
    ax.set_xlabel('Number of Pharmacogene Connections', fontsize=12, fontweight='bold')
    ax.set_title('MHT Compound Pharmacogene Connectivity in SPOKE', fontsize=14, fontweight='bold')
    for i, val in enumerate(pivot['total'].values):
        ax.text(val + 0.1, i, str(int(val)), va='center', fontsize=9, fontweight='bold')
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(True, alpha=0.3, axis='x')
    sns.despine(left=True)
    plt.tight_layout()
    if save_path:
        fig.savefig(f'{save_path}pharmacogene_mht_connectivity.pdf', dpi=300, bbox_inches='tight')
        fig.savefig(f'{save_path}pharmacogene_mht_connectivity.png', dpi=300, bbox_inches='tight')
        print(f'✓ Connectivity chart saved to {save_path}')
    return fig

def run_pharmacogene_analysis(open_spoke, SPOKE_User, SPOKE_Password, SPOKE_DEV_URI, save_path=''):
    """
    Main entry point.

    Usage:
        all_edges, matrix = run_pharmacogene_analysis(
            open_spoke, SPOKE_User, SPOKE_Password, SPOKE_DEV_URI,
            save_path='${SPOKE_DATA_ROOT}/mht_results/'
        )
        plt.show()
    """
    all_edges = query_spoke_pharmacogene_mht(open_spoke, SPOKE_User, SPOKE_Password, SPOKE_DEV_URI)
    matrix = build_interaction_matrix(all_edges)
    fig1 = plot_heatmap(matrix, save_path)
    fig2, gene_table_df = plot_gene_table(all_edges, save_path)
    fig3 = plot_compound_connectivity(all_edges, save_path)
    all_edges.to_csv(f'{save_path}pharmacogene_mht_edges.csv', index=False)
    matrix.to_csv(f'{save_path}pharmacogene_mht_matrix.csv')
    gene_table_df.to_csv(f'{save_path}pharmacogene_table.csv', index=False)
    print(f'\n✓ All results saved to {save_path}')
    return (all_edges, matrix)
if __name__ == '__main__':
    print('\n    ================================================================\n    Pharmacogene × MHT Compound Interaction Analysis\n    ================================================================\n\n    Outputs:\n        1. pharmacogene_mht_heatmap.png/pdf      — Annotated heatmap\n        2. pharmacogene_table.png/pdf             — Gene function table\n        3. pharmacogene_mht_connectivity.png/pdf  — Bar chart\n        4. pharmacogene_mht_edges.csv             — Raw edge table\n        5. pharmacogene_mht_matrix.csv            — Interaction matrix\n        6. pharmacogene_table.csv                 — Gene table data\n    ')
