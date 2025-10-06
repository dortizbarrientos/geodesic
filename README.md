## GitHub Repository Structure

```
senecio-gmatrix-affine-analysis/
├── README.md
├── LICENSE
├── requirements.txt
├── data/
│   ├── README.md
│   ├── senecio_bootstrap_results.csv
│   └── senecio_geodesic_path.csv
├── scripts/
│   ├── senecio_complete_analysis.py
│   ├── senecio_affine_analysis_2.py
│   └── senecio_figure1d_clean.py
├── figures/
│   └── (output figures will be saved here)
└── docs/
    ├── METHODS.md
    └── INTERPRETATION.md
```

## Main README.md Content

# Affine-Invariant Analysis of G-Matrix Evolution in Senecio lautus

Geometric analysis of genetic variance-covariance matrix evolution across four ecotypes of *Senecio lautus*demonstrates how adaptation to contrasting coastal environments reorganizes the G-matrix structure.

## Overview

This repository contains Python implementations of affine-invariant distance methods applied to G-matrices from Walter et al. (2018, *American Naturalist*). The analysis quantifies G-matrix divergence by decomposing total distance into three orthogonal components:

- **Scale**: Changes in overall genetic variance magnitude
- **Anisotropy**: Changes in integration strength (variance distribution across eigenvalues)
- **Orientation**: Changes in eigenvector directions (trait correlations)

## Key Findings

1. **Headland's unique architecture**: The prostrate headland ecotype shows 1.7× higher anisotropy distances compared to upright ecotypes, driven by extreme integration (50% variance in *g*max vs. 24-32% in other ecotypes)

2. **Orientation dominates divergence**: Across all pairwise comparisons, orientation changes account for 76.5% of total distance, but anisotropy contributes substantially (11.1%)—unlike other systems where anisotropy is negligible

3. **Geodesic analysis**: The shortest path from Dune to Headland shows a smooth, monotonic increase in integration strength, suggesting coordinated evolutionary reorganization

## Scripts

### Primary Analysis: `senecio_complete_analysis.py`

Complete workflow including distance calculations, bootstrap confidence intervals, geodesic path analysis, and publication-quality visualizations.

**Outputs:**
- 6-panel figure with distance matrices, component decomposition, and geodesic analysis
- `senecio_bootstrap_results.csv`: Bootstrap CIs for all pairwise distances
- `senecio_geodesic_path.csv`: Properties along Dune→Headland geodesic

**Runtime:** ~2-3 minutes for 200 bootstrap replicates

### Simplified Analysis: `senecio_affine_analysis_2.py`

Streamlined script for calculating affine-invariant distances without bootstrap or geodesics. Use this for quick calculations or to understand the core method.

### Figure Generation: `senecio_figure1d_clean.py`

Generates Figure 1D showing *g*max trait composition across ecotypes. Minimal dependencies, optimized for PNAS style.

## Installation

```bash
# Clone repository
git clone https://github.com/[your-username]/senecio-gmatrix-affine-analysis.git
cd senecio-gmatrix-affine-analysis

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Requirements

```
numpy>=1.21.0
pandas>=1.3.0
matplotlib>=3.4.0
seaborn>=0.11.0
scipy>=1.7.0
scikit-learn>=0.24.0
```

## Usage

### Run complete analysis:

```python
python scripts/senecio_complete_analysis.py
```

### Run simplified distance calculations:

```python
python scripts/senecio_affine_analysis_2.py
```

### Generate Figure 1D only:

```python
python scripts/senecio_figure1d_clean.py
```

## Data Source

G-matrices extracted from Walter et al. (2018) Table S1. The G-matrices are hard-coded in the scripts as they represent published posterior means from Bayesian half-sib analyses.

**Reference:**
Walter, G. M., Aguirre, J. D., Blows, M. W., & Ortiz-Barrientos, D. (2018). Evolution of genetic variance during adaptive radiation. *The American Naturalist*, 191(4), E108-E128.

## Method Details

### Affine-Invariant Distance

The affine-invariant (Riemannian) distance between G-matrices **G₁** and **G₂**:

```
d(G₁, G₂) = ||log(G₁^(-1/2) G₂ G₁^(-1/2))||_F
```

This metric is invariant to:
- Linear transformations of trait space
- Changes in measurement units
- Trait rotations

### Component Decomposition

Total distance decomposes via the Pythagorean theorem on the manifold:

```
d²_total = d²_scale + d²_aniso + d²_orient
```

**Scale component:**
```
d_scale = √n |log(s₂/s₁)|
```
where s = geometric mean eigenvalue

**Anisotropy component:**
```
d_aniso = √Σ[log(λ₂ᵢ/s₂) - log(λ₁ᵢ/s₁)]²
```
where λᵢ are eigenvalues

**Orientation component:** Residual distance after accounting for scale and anisotropy

### Bootstrap Procedure

Parametric bootstrap (n=1000 replicates):
1. Add proportional noise to eigenvalues (15% CV)
2. Reconstruct G-matrix via eigendecomposition
3. Calculate distances for each replicate
4. Extract 95% confidence intervals (2.5th, 97.5th percentiles)

### Geodesic Path

The geodesic from **G₀** to **G₁** parameterized by t ∈ [0,1]:

```
G(t) = G₀^(1/2) [G₀^(-1/2) G₁ G₀^(-1/2)]^t G₀^(1/2)
```

Computed via eigendecomposition of the matrix logarithm.

## Comparison to Tensor Analysis

Walter et al. (2018) used covariance tensor analysis (Aguirre et al. 2014), which identified three eigentensors explaining G-matrix divergence. The first eigentensor (E₁) explained 38% of variance, dominated by architectural trait orientations.

The affine-invariant framework provides complementary information:
- **Agreement**: Both methods identify orientation changes as primary driver; Headland diverges most
- **Advantage**: Affine decomposition quantifies integration strength changes (anisotropy) specifically, which tensor analysis does not separate explicitly
- **Result**: Anisotropy accounts for 11.1% of total distance—detectable unlike in other systems (e.g., *Anolis* where <1%)

## Citation

If you use this code, please cite:

1. **Original data source:**  
   Walter, G. M., Aguirre, J. D., Blows, M. W., & Ortiz-Barrientos, D. (2018). Evolution of genetic variance during adaptive radiation. *The American Naturalist*, 191(4), E108-E128.

2. **Affine-invariant methods:**  
   Ortiz-Barrientos, D., Burrage, K., and Burrage, P. (2025). Adaptive reorganization of the G-matrix on a curved space. *Submitted*, BioRxiv: TBA

## License

MIT License - see LICENSE file

## Contact

[Your name and email]

## Acknowledgments

Analysis builds on theoretical framework developed by [cite relevant geometric statistics papers]. Thanks to Greg Walter for making G-matrix data publicly available.
```

## requirements.txt

```
numpy==1.24.3
pandas==2.0.3
matplotlib==3.7.2
seaborn==0.12.2
scipy==1.11.1
scikit-learn==1.3.0
```

## docs/METHODS.md

Create a detailed mathematical explanation of:
1. Why affine-invariance matters for G-matrices
2. Geometric interpretation of the paraboloid framework
3. Proof that components are orthogonal
4. Comparison to other distance metrics (Euclidean, Procrustes, tensor)

## docs/INTERPRETATION.md

Biological interpretation guide:
1. What does anisotropy mean evolutionarily?
2. When should we expect anisotropy to change vs. stay constant?
3. How to interpret geodesic paths
4. Limitations and assumptions

## Additional Files Needed

1. **LICENSE** - Choose MIT, GPL, or BSD depending on your preference
2. **CHANGELOG.md** - Track version updates
3. **.gitignore** - Standard Python gitignore

This structure provides clear navigation, reproducible workflows, and sufficient documentation for both users who want to run the analysis and those who want to understand the methods deeply.
