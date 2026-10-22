# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def plot_standalone_heatmap_legend(save_path='', filename='pharmacogene_mht_heatmap_legend', include_all_three=True):
    """
    Standalone legend figure for the Pharmacogene × MHT heatmap.
    Saves PNG and PDF if save_path is provided.
    """
    edge_type_colors = {'Binds protein': '#3B82F6', 'Interacts with protein': '#10B981', 'Binds + Interacts': '#8B5CF6', 'Upregulates gene': '#F97316', 'Downregulates gene': '#EF4444', 'Binding + Regulation': '#1E3A5F'}
    if include_all_three:
        edge_type_colors['Binds + Interacts + Regulation'] = '#111827'
    compound_class_colors = {'Estrogen': '#E74C8B', 'Estrogen (CE)': '#F4A6C9', 'Progestin': '#6366F1', 'SERM': '#F59E0B'}
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.axis('off')
    ax.text(0.5, 0.95, 'Pharmacogene × MHT Interaction Legend', ha='center', va='top', fontsize=20, fontweight='bold')
    ax.text(0.08, 0.8, 'Edge Type', ha='left', va='center', fontsize=16, fontweight='bold')
    ax.text(0.58, 0.8, 'Drug Class', ha='left', va='center', fontsize=16, fontweight='bold')
    y_start = 0.68
    y_step = 0.095
    for i, (label, color) in enumerate(edge_type_colors.items()):
        y = y_start - i * y_step
        rect = mpatches.FancyBboxPatch((0.08, y - 0.025), 0.06, 0.05, boxstyle='round,pad=0.004,rounding_size=0.008', facecolor=color, edgecolor='black', linewidth=0.6, transform=ax.transAxes)
        ax.add_patch(rect)
        ax.text(0.16, y, label, ha='left', va='center', fontsize=13.5, fontweight='bold')
    y_start_drug = 0.68
    for i, (label, color) in enumerate(compound_class_colors.items()):
        y = y_start_drug - i * y_step
        rect = mpatches.FancyBboxPatch((0.58, y - 0.025), 0.06, 0.05, boxstyle='round,pad=0.004,rounding_size=0.008', facecolor=color, edgecolor='black', linewidth=0.6, transform=ax.transAxes)
        ax.add_patch(rect)
        ax.text(0.66, y, label, ha='left', va='center', fontsize=13.5, fontweight='bold')
    plt.tight_layout()
    if save_path:
        fig.savefig(f'{save_path}{filename}.png', dpi=300, bbox_inches='tight', transparent=True)
        fig.savefig(f'{save_path}{filename}.pdf', dpi=300, bbox_inches='tight', transparent=True)
        print(f'✓ Legend saved to: {save_path}{filename}.png/pdf')
    return fig
fig = plot_standalone_heatmap_legend(save_path='${SPOKE_DATA_ROOT}/mht_results/')
plt.show()
