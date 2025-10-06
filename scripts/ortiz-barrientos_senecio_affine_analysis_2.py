"""
Senecio G-Matrix Affine-Invariant Analysis - CORRECTED VERSION
Complete implementation for Walter et al. (2018) data

Author: Analysis based on Walter et al. 2018 Am Nat
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.linalg import sqrtm, logm
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CORRECTED: DIRECT EXTRACTION FROM TABLE S1 COVARIANCE VALUES
# ============================================================================

# The table provides variances on diagonal and COVARIANCES in upper triangle
# We'll extract directly from the published covariance matrices

def create_g_from_published_table(upper_triangle_covs, diagonal_vars):
    """
    Create G-matrix directly from published covariances.
    Table S1 shows variances on diagonal, covariances above diagonal.
    """
    n = len(diagonal_vars)
    G = np.zeros((n, n))
    
    # Set diagonal
    np.fill_diagonal(G, diagonal_vars)
    
    # Set upper triangle
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            G[i, j] = upper_triangle_covs[idx]
            G[j, i] = upper_triangle_covs[idx]  # Symmetric
            idx += 1
    
    return G

# Extract from Table S1 - reading across rows for upper triangle
# Dune: variances on diagonal, covariances above diagonal
dune_vars = [0.294, 0.196, 0.366, 0.399, 0.058, 0.126, 0.060, 0.203, 0.197, 0.088]
dune_covs = [
    # Row 1: covs (1,2) through (1,10)
    0.096, 0.060, 0.067, 0.013, -0.019, 0.007, 0.019, -0.021, -0.012,
    # Row 2: covs (2,3) through (2,10)  
    0.049, 0.011, 0.013, 0.009, 0.001, 0.009, -0.016, 0.005,
    # Row 3: covs (3,4) through (3,10)
    0.005, 0.001, -0.056, 0.014, 0.014, 0.020, -0.024,
    # Row 4: covs (4,5) through (4,10)
    0.001, -0.009, 0.013, 0.003, 0.008, -0.005,
    # Row 5: covs (5,6) through (5,10)
    -0.007, 0.005, 0.019, -0.023, 0.005,
    # Row 6: covs (6,7) through (6,10)
    -0.009, -0.026, -0.010, 0.040,
    # Row 7: covs (7,8) through (7,10)
    0.019, -0.019, -0.008,
    # Row 8: covs (8,9) through (8,10)
    -0.126, -0.067,
    # Row 9: covs (9,10)
    0.027
]

# Headland
headland_vars = [0.421, 0.737, 0.768, 0.822, 1.003, 0.071, 0.083, 0.219, 0.158, 0.176]
headland_covs = [
    0.155, -0.110, 0.286, 0.294, -0.054, 0.062, -0.121, -0.015, 0.122,
    0.171, 0.398, 0.493, -0.064, 0.062, -0.181, 0.013, 0.178,
    -0.068, -0.127, -0.003, -0.021, 0.092, -0.027, -0.047,
    0.585, -0.059, 0.107, -0.242, -0.018, 0.233,
    -0.097, 0.112, -0.310, 0.019, 0.276,
    -0.020, 0.020, -0.009, -0.010,
    -0.024, -0.019, 0.038,
    -0.063, -0.123,
    -0.003
]

# Tableland  
tableland_vars = [0.342, 0.455, 0.242, 0.812, 0.236, 0.086, 0.266, 0.165, 0.262, 0.514]
tableland_covs = [
    0.193, 0.028, 0.243, 0.022, -0.035, 0.019, 0.007, 0.001, 0.039,
    0.013, 0.203, 0.102, -0.037, 0.108, -0.038, -0.023, 0.162,
    0.024, -0.020, 0.004, -0.026, 0.001, 0.026, -0.082,
    -0.014, -0.034, -0.018, 0.062, -0.024, -0.032,
    -0.033, 0.122, -0.047, -0.099, 0.188,
    -0.024, -0.012, -0.003, -0.006,
    -0.041, -0.065, 0.202,
    -0.057, -0.142,
    -0.084
]

# Woodland
woodland_vars = [0.173, 0.106, 0.278, 0.135, 0.084, 0.162, 0.095, 0.269, 0.151, 0.094]
woodland_covs = [
    0.029, -0.011, 0.015, 0.014, 0.000, 0.001, -0.027, 0.017, 0.000,
    0.015, 0.017, -0.009, 0.024, -0.006, -0.025, 0.036, 0.004,
    0.008, -0.031, -0.027, -0.009, 0.061, -0.010, -0.026,
    0.017, -0.037, 0.017, 0.022, -0.019, -0.022,
    -0.028, 0.037, -0.013, -0.022, -0.012,
    -0.016, -0.072, 0.035, 0.033,
    -0.034, -0.048, -0.018,
    -0.053, -0.045,
    0.034
]

# Construct G-matrices
G_dune = create_g_from_published_table(dune_covs, dune_vars)
G_headland = create_g_from_published_table(headland_covs, headland_vars)
G_tableland = create_g_from_published_table(tableland_covs, tableland_vars)
G_woodland = create_g_from_published_table(woodland_covs, woodland_vars)

G_matrices = {
    'Dune': G_dune,
    'Headland': G_headland,
    'Tableland': G_tableland,
    'Woodland': G_woodland
}

trait_names = ['Height', 'StemL/W', 'Branches', 'StemDiam', 'LeafArea', 
               'Perim²/Area²', 'Circularity', 'Nindents', 'IndentWidth', 'IndentDepth']

# ============================================================================
# UTILITY FUNCTIONS - CORRECTED
# ============================================================================

def affine_invariant_distance(G1, G2, regularize=True, epsilon=1e-8):
    """
    Calculate affine-invariant distance with optional regularization.
    """
    if regularize:
        # Add small regularization to ensure positive definiteness
        G1 = G1 + epsilon * np.eye(len(G1))
        G2 = G2 + epsilon * np.eye(len(G2))
    
    try:
        G1_sqrt = sqrtm(G1).real
        G1_inv_sqrt = np.linalg.inv(G1_sqrt)
        M = G1_inv_sqrt @ G2 @ G1_inv_sqrt
        log_M = logm(M).real
        distance = np.linalg.norm(log_M, 'fro')
        return distance
    except:
        return np.nan

def decompose_distance(G1, G2, regularize=True):
    """
    Decompose affine distance into components with error handling.
    """
    epsilon = 1e-8
    if regularize:
        G1 = G1 + epsilon * np.eye(len(G1))
        G2 = G2 + epsilon * np.eye(len(G2))
    
    try:
        d_total = affine_invariant_distance(G1, G2, regularize=False)
        
        evals1, evecs1 = np.linalg.eigh(G1)
        evals2, evecs2 = np.linalg.eigh(G2)
        
        # Ensure positive eigenvalues
        evals1 = np.maximum(evals1, epsilon)
        evals2 = np.maximum(evals2, epsilon)
        
        # Scale component
        scale1 = np.exp(np.mean(np.log(evals1)))
        scale2 = np.exp(np.mean(np.log(evals2)))
        n = len(evals1)
        d_scale = np.sqrt(n) * abs(np.log(scale2) - np.log(scale1))
        
        # Anisotropy component
        norm_evals1 = evals1 / scale1
        norm_evals2 = evals2 / scale2
        d_aniso = np.sqrt(np.sum((np.log(norm_evals2) - np.log(norm_evals1))**2))
        
        # Orientation component
        d_orient_sq = max(0, d_total**2 - d_scale**2 - d_aniso**2)
        d_orient = np.sqrt(d_orient_sq)
        
        total_sq = d_total**2
        return {
            'total': d_total,
            'scale': d_scale,
            'anisotropy': d_aniso,
            'orientation': d_orient,
            'scale_pct': (d_scale**2 / total_sq) * 100 if total_sq > 0 else 0,
            'aniso_pct': (d_aniso**2 / total_sq) * 100 if total_sq > 0 else 0,
            'orient_pct': (d_orient**2 / total_sq) * 100 if total_sq > 0 else 0
        }
    except Exception as e:
        print(f"Error in decomposition: {e}")
        return {
            'total': np.nan, 'scale': np.nan, 'anisotropy': np.nan, 'orientation': np.nan,
            'scale_pct': np.nan, 'aniso_pct': np.nan, 'orient_pct': np.nan
        }

# ============================================================================
# VALIDATION
# ============================================================================

print("="*80)
print("SENECIO G-MATRIX AFFINE-INVARIANT ANALYSIS")
print("Data from Walter et al. (2018) American Naturalist")
print("="*80)

print("\nStep 1: VALIDATION")
print("-"*80)

epsilon = 1e-8
for name, G in G_matrices.items():
    # Check properties
    is_sym = np.allclose(G, G.T)
    evals = np.linalg.eigvalsh(G)
    
    # Check for small/negative eigenvalues
    min_eval = np.min(evals)
    is_pd = np.all(evals > 0)
    
    # Calculate % variance in gmax
    total_var = np.sum(evals[evals > 0])  # Only positive eigenvalues
    pct_gmax = (evals[-1] / total_var) * 100
    
    print(f"\n{name} ecotype:")
    print(f"  Symmetric: {is_sym}")
    print(f"  Min eigenvalue: {min_eval:.6f}")
    print(f"  All positive: {is_pd}")
    if not is_pd:
        print(f"  → Will add regularization (ε={epsilon})")
    print(f"  % variance in gmax: {pct_gmax:.1f}%")
    print(f"  Trace (total variance): {np.trace(G):.3f}")
    print(f"  Top 3 eigenvalues: {evals[-3:][::-1]}")

print("\n" + "-"*80)
print("Expected from Table 1:")
print("  Dune=24%, Headland=50%, Tableland=32%, Woodland=26%")

# ============================================================================
# PAIRWISE DISTANCES
# ============================================================================

print("\n" + "="*80)
print("Step 2: PAIRWISE AFFINE-INVARIANT DISTANCES")
print("="*80)

ecotypes = list(G_matrices.keys())
results = []

print("\nCalculating pairwise distances with regularization...")

for eco1, eco2 in combinations(ecotypes, 2):
    G1 = G_matrices[eco1]
    G2 = G_matrices[eco2]
    
    decomp = decompose_distance(G1, G2, regularize=True)
    
    results.append({
        'Ecotype 1': eco1,
        'Ecotype 2': eco2,
        'Total': decomp['total'],
        'Scale': decomp['scale'],
        'Anisotropy': decomp['anisotropy'],
        'Orientation': decomp['orientation'],
        'Scale %': decomp['scale_pct'],
        'Aniso %': decomp['aniso_pct'],
        'Orient %': decomp['orient_pct']
    })

df = pd.DataFrame(results)

print("\nTotal distances:")
print("-"*80)
for _, row in df.iterrows():
    print(f"  {row['Ecotype 1']:12s} vs {row['Ecotype 2']:12s}: {row['Total']:.4f}")

print("\n\nComponent values:")
print("-"*80)
print(f"{'Comparison':<28} {'Total':>8} {'Scale':>8} {'Aniso':>8} {'Orient':>8}")
print("-"*80)
for _, row in df.iterrows():
    comp = f"{row['Ecotype 1']} vs {row['Ecotype 2']}"
    print(f"{comp:<28} {row['Total']:8.4f} {row['Scale']:8.4f} "
          f"{row['Anisotropy']:8.4f} {row['Orientation']:8.4f}")

print("\n\nComponent percentages:")
print("-"*80)
print(f"{'Comparison':<28} {'Scale':>10} {'Aniso':>10} {'Orient':>10}")
print("-"*80)
for _, row in df.iterrows():
    comp = f"{row['Ecotype 1']} vs {row['Ecotype 2']}"
    print(f"{comp:<28} {row['Scale %']:9.1f}% {row['Aniso %']:9.1f}% {row['Orient %']:9.1f}%")

# ============================================================================
# HYPOTHESIS TESTING
# ============================================================================

print("\n" + "="*80)
print("Step 3: BIOLOGICAL HYPOTHESES")
print("="*80)

# H1: Headland comparisons
headland_rows = df[(df['Ecotype 1'] == 'Headland') | (df['Ecotype 2'] == 'Headland')]

print("\nH1: Headland divergence patterns")
print("-"*80)
for _, row in headland_rows.iterrows():
    other = row['Ecotype 1'] if row['Ecotype 2'] == 'Headland' else row['Ecotype 2']
    print(f"  Headland vs {other:12s}: Aniso={row['Aniso %']:5.1f}%, Orient={row['Orient %']:5.1f}%")

# H2: Upright comparisons
upright_rows = df[df.apply(
    lambda r: all(e in ['Dune', 'Tableland', 'Woodland'] 
                  for e in [r['Ecotype 1'], r['Ecotype 2']]), axis=1)]

print("\nH2: Within-upright comparisons")
print("-"*80)
for _, row in upright_rows.iterrows():
    print(f"  {row['Ecotype 1']:12s} vs {row['Ecotype 2']:12s}: "
          f"Aniso={row['Aniso %']:5.1f}%, Orient={row['Orient %']:5.1f}%")

# Summary stats
print("\nSummary statistics:")
print("-"*80)
print(f"Mean contribution across all pairs:")
print(f"  Orientation: {df['Orient %'].mean():.1f}%")
print(f"  Anisotropy:  {df['Aniso %'].mean():.1f}%")
print(f"  Scale:       {df['Scale %'].mean():.1f}%")

if not headland_rows['Anisotropy'].isna().all() and not upright_rows['Anisotropy'].isna().all():
    mean_h = headland_rows['Anisotropy'].mean()
    mean_u = upright_rows['Anisotropy'].mean()
    print(f"\nAnisotropy distances:")
    print(f"  Headland mean: {mean_h:.4f}")
    print(f"  Upright mean:  {mean_u:.4f}")
    print(f"  Ratio: {mean_h/mean_u:.2f}x")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
print("\nNote: If % gmax values still don't match Table 1, there may be")
print("      discrepancies in how the published values were calculated.")
print("      Proceed with comparing PATTERNS rather than exact values.")
