"""
Figure 1D: g_max Eigenvector Loadings - PNAS Clean Style
==========================================================

Simplified, clean version optimized for PNAS submission.
Uses minimal colors and clear visual hierarchy.

Author: Daniel Ortiz-Barrientos
Date: October 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ============================================================================
# DATA
# ============================================================================

trait_names = [
    'Height', 'Stem L/W', 'Branches', 'Stem Diam', 'Leaf Area',
    'Perim²/Area²', 'Circularity', 'N Indents', 'Indent W', 'Indent D'
]

loadings = {
    'Dune': [0.35, 0.21, 0.59, 0.66, 0.03, -0.06, 0.07, 0.15, -0.15, -0.01],
    'Headland': [0.26, 0.42, -0.07, 0.52, 0.61, -0.06, 0.09, -0.23, 0.00, 0.21],
    'Tableland': [0.39, 0.48, 0.01, 0.71, 0.15, -0.06, 0.15, 0.01, -0.09, 0.24],
    'Woodland': [-0.15, -0.11, 0.51, 0.04, -0.02, -0.35, -0.03, 0.68, -0.27, -0.22]
}

pct_gmax = {'Dune': 24.4, 'Headland': 49.7, 'Tableland': 31.7, 'Woodland': 26.0}

# ============================================================================
# PLOTTING
# ============================================================================

def plot_panel_clean(ax, ecotype, values, pct, is_headland=False):
    """Simple, clean bar plot for one ecotype."""
    
    n = len(values)
    y_pos = np.arange(n)
    
    # Simple two-color scheme: positive (blue), negative (red)
    # Headland uses red theme, others use blue
    colors = []
    for val in values:
        if is_headland:
            colors.append('#e74c3c' if val > 0 else '#c0392b')
        else:
            colors.append('#3498db' if val > 0 else '#2980b9')
    
    # Create bars with emphasis on strong loadings (|val| > 0.25)
    bars = ax.barh(y_pos, values, height=0.7, color=colors, 
                   edgecolor='black', linewidth=0.5, alpha=0.9)
    
    # Emphasize strong loadings with thicker edge
    for i, (bar, val) in enumerate(zip(bars, values)):
        if abs(val) > 0.25:
            bar.set_linewidth(1.5)
            bar.set_alpha(1.0)
    
    # Zero line
    ax.axvline(0, color='black', linewidth=1, zorder=0)
    
    # Formatting
    ax.set_yticks(y_pos)
    ax.set_yticklabels(trait_names, fontsize=8)
    ax.set_xlim(-0.75, 0.75)
    ax.set_xticks([-0.6, -0.3, 0, 0.3, 0.6])
    ax.tick_params(axis='x', labelsize=8)
    ax.grid(axis='x', alpha=0.2, linestyle='--')
    ax.set_axisbelow(True)
    ax.invert_yaxis()
    
    # Title
    color = '#e74c3c' if is_headland else '#2c3e50'
    ax.set_title(f'{ecotype} ({pct:.1f}%)', fontsize=10, 
                fontweight='bold', color=color, pad=8)
    
    # Only show x-label on bottom panels
    if ecotype in ['Tableland', 'Woodland']:
        ax.set_xlabel('Loading', fontsize=9, fontweight='bold')

# ============================================================================
# CREATE FIGURE
# ============================================================================

fig = plt.figure(figsize=(10, 8))
gs = GridSpec(2, 2, hspace=0.3, wspace=0.25, 
              left=0.12, right=0.96, top=0.93, bottom=0.1)

# Create panels
ecotypes = ['Dune', 'Headland', 'Tableland', 'Woodland']
positions = [(0,0), (0,1), (1,0), (1,1)]

for eco, (row, col) in zip(ecotypes, positions):
    ax = fig.add_subplot(gs[row, col])
    is_headland = (eco == 'Headland')
    plot_panel_clean(ax, eco, loadings[eco], pct_gmax[eco], is_headland)

# Overall title
fig.text(0.5, 0.97, 'D) g_max Trait Composition', 
         ha='center', fontsize=12, fontweight='bold')

# ============================================================================
# SAVE
# ============================================================================

plt.savefig('Figure_1D_clean.png', dpi=300, bbox_inches='tight', 
            facecolor='white')
plt.savefig('Figure_1D_clean.pdf', bbox_inches='tight')

print("✓ Clean version saved!")
print("\nFigure features:")
print("  • Minimal color scheme (blue for upright, red for headland)")
print("  • Strong loadings (|x| > 0.25) have thicker borders")
print("  • Optimized for PNAS multi-panel figure")
print("\nBiological pattern:")
print("  Headland: Coordinated size increase (Leaf Area 0.61, Stem Diam 0.52)")
print("  Upright: Variable g_max reflecting different selective environments")

plt.show()
