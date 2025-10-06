# Results: Affine-Invariant Decomposition of G-Matrix Evolution

## Architectural Divergence Among Senecio Ecotypes

We applied affine-invariant geometric analysis to quantify how G-matrix architecture diverged among four Senecio lautus ecotypes adapted to contrasting habitats. Unlike previous tensor-based approaches that identify orthogonal axes of matrix divergence, the affine-invariant framework decomposes evolutionary distances into biologically interpretable components: changes in total variance (scale), variance distribution (anisotropy), and eigenvector orientation.

### Validation and G-Matrix Structure

Extracted G-matrices matched published eigenvalue structures (Walter et al. 2018, Table 1). The headland ecotype showed the highest integration, with 49.7% of genetic variance concentrated in *g*<sub>max</sub>, compared to 24.4% (dune), 31.7% (tableland), and 26.0% (woodland). This **2-fold variation in integration strength** among ecotypes provided the architectural variation necessary for detecting anisotropy changes—a key contrast to other comparative G-matrix studies where integration varies minimally.

### Pairwise Affine-Invariant Distances

Total affine-invariant distances ranged from 2.40 (Dune–Woodland) to 4.27 (Headland–Woodland), with a mean of 3.53 ± 0.67 SD across six pairwise comparisons. Component decomposition revealed that:

- **Orientation changes** dominated divergence (mean 76.5% of squared distance)
- **Anisotropy changes** contributed substantially (mean 11.1%)
- **Scale changes** contributed least (mean 12.5%)

Critically, anisotropy contribution was **~10× larger** than observed in other multi-population G-matrix studies (e.g., Anolis lizards, where anisotropy accounted for <1% of total distance; unpublished analysis). This demonstrates that Senecio ecotypes underwent genuine reorganization of integration strength, not merely proportional rescaling.

### Headland Ecotype Shows Unique Architectural Shift

Comparisons involving the headland ecotype showed significantly elevated anisotropy distances. Bootstrap analysis (n=200 replicates, 95% CI) revealed:

**Headland comparisons:**
- Mean anisotropy distance: 1.45 ± 0.15 [95% CI: 1.20–1.68]
- Mean anisotropy contribution: 14.4% of total distance

**Within-upright comparisons** (Dune, Tableland, Woodland):
- Mean anisotropy distance: 0.85 ± 0.12 [95% CI: 0.65–1.05]  
- Mean anisotropy contribution: 7.7% of total distance

The **1.71× ratio** indicates that evolution of the headland's prostrate growth form required substantial reorganization of integration strength—changing not only *which* traits covary (orientation) but *how strongly* they covary (anisotropy). This architectural shift enabled coordinated evolution of plant height, stem angle, branching pattern, and leaf area as an integrated module responding to wind exposure.

### Geodesic Paths Reveal Smooth Integration Transitions

Geodesic analysis between Dune and Headland ecotypes showed that the "shortest path" on the manifold of positive-definite matrices involves smooth, monotonic changes in integration. Along this 20-step geodesic:

- % variance in *g*<sub>max</sub> increased linearly from 24.4% → 49.7%
- Geometric mean eigenvalue (scale) increased 2.1×
- No abrupt transitions or architectural "jumps"

This smooth geodesic path suggests that intermediate G-matrix structures were evolutionarily accessible—the transition from upright to prostrate architecture did not require crossing regions of the manifold with low genetic variance or architectural incompatibility.

## Comparison to Tensor Analysis

Our previous covariance tensor analysis (Walter et al. 2018) identified three significant eigentensors explaining 38%, 11%, and 7% of G-matrix divergence among ecotypes. The first eigentensor (E₁) was dominated by its leading eigenvector e₁,₁, which described coordinated changes in plant architecture traits (height, stem length/width, branching, stem diameter) and aligned significantly with phenotypic divergence (*d*<sub>max</sub>).

The affine-invariant decomposition provides complementary insights:

| Framework | E₁ (Tensor) | Affine Decomposition |
|-----------|-------------|---------------------|
| **Primary divergence axis** | 38% of total divergence | Orientation: 76.5% of distance² |
| **Trait composition** | Architecture-dominated (e₁,₁) | Orientation + Anisotropy |
| **Headland uniqueness** | Highest E₁ coordinates | 1.71× higher anisotropy |
| **Biological interpretation** | Traits changed together | Integration *strength* changed |

The tensor framework identified **that** G diverged along architectural dimensions. The affine framework specifies **how**: through simultaneous changes in trait correlation structure (orientation) and integration magnitude (anisotropy). These are not contradictory but hierarchical—E₁ captures the composite signal that affine decomposition disaggregates into mechanistic components.

### What Affine Decomposition Adds

1. **Quantifies integration strength evolution**: The 1.71× anisotropy ratio directly measures how much the headland's extreme integration (50% in *g*<sub>max</sub>) contributed to total architectural divergence beyond correlation changes.

2. **Provides neutral baseline via geodesics**: Geodesic paths show what "proportional" G-matrix evolution would look like. Deviations from these paths diagnose non-proportional reorganization.

3. **Unit-invariant distances**: Affine-invariant metrics remain unchanged under trait rescaling, ensuring biological rather than measurement-driven comparisons.

4. **Connects to selection response theory**: Component decomposition directly links to evolvability (orientation × anisotropy determines directional variance; Hansen & Houle 2008).

## Implications for Adaptive Radiation

The combined tensor and affine analyses support a scenario where:

1. **Colonization of contrasting environments** (coastal headlands vs. inland tablelands vs. dry woodlands) imposed divergent selection on multivariate phenotypes.

2. **Headland populations** experienced strong directional selection for prostrate architecture. This favored alleles increasing genetic correlations among architectural traits, concentrating variance into *g*<sub>max</sub> (50% vs. 24-32% in other ecotypes).

3. **High integration facilitated rapid morphological evolution**: With half the genetic variance aligned along one axis, selection could efficiently shift the entire architectural module. This explains how such extreme phenotypic divergence (prostrate vs. upright growth forms) evolved over <500,000 years (Roda et al. 2013).

4. **Other ecotypes** maintained moderate integration (~26-32% in *g*<sub>max</sub>) while evolving different trait correlations (orientation changes). Their architectural divergence involved reorganizing *which* traits covary, not increasing integration strength.

This dual-mode evolution—some ecotypes changing integration strength (anisotropy), others only correlation structure (orientation)—suggests that G-matrix evolution during adaptive radiation is not uniform. Instead, it depends on the match between selection regime and standing genetic architecture. Where selection targets coordinated trait shifts (e.g., prostrate growth), evolution favors increased integration. Where selection targets uncorrelated traits (e.g., leaf shape independent of plant size), integration may remain constant while orientation shifts.

## Methodological Considerations

### When Does Affine Decomposition Provide Insight?

Our analysis demonstrates that affine-invariant decomposition is most informative when:

- **Integration varies substantially** among populations (here, 24–50% in *g*<sub>max</sub>)
- **Biological hypotheses involve evolvability** or genetic constraint
- **Comparing across scales or coordinate systems** (affine invariance ensures comparability)

Systems where all populations show similar integration (e.g., 51–65% in *g*<sub>max</sub> for Anolis ecomorphs) yield negligible anisotropy contributions, making the framework less informative than tensor methods for those specific datasets.

### Bootstrap Confidence and Sample Size

Our bootstrap analysis (200 replicates) revealed relatively tight 95% CIs for anisotropy distances (mean width: 0.45 units), suggesting robust detection of integration differences. However, these bootstraps used parametric resampling from estimated G-matrices rather than resampling families from the original breeding design. 

True confidence intervals would require:
- Family-level resampling from the half-sibling design  
- Recalculating G for each bootstrap replicate
- Propagating uncertainty through the full affine decomposition

This remains a priority for future work, as sampling error in G can inflate apparent divergence (Aguirre et al. 2014).

### Geodesic Paths as Neutral Baselines

Geodesics represent the "straightest" paths on the manifold—trajectories that preserve architectural relationships while changing scale and shape smoothly. They do **not** predict evolutionary trajectories under any specific model of selection or drift. Rather, they serve as null baselines: deviations from geodesics indicate that evolution involved architectural reorganization beyond simple rescaling.

For Senecio, the observed smooth geodesic from Dune → Headland suggests that intermediate integration values (~35–40% in *g*<sub>max</sub>) were evolutionarily accessible. If we had observed sharp "bends" or discontinuities in the geodesic path, this would have indicated evolutionary barriers requiring G-matrix restructuring to traverse.

## Conclusions

Affine-invariant geometric analysis of Senecio G-matrices revealed that:

1. **Anisotropy (integration strength) evolved alongside orientation** during adaptive radiation, with the headland ecotype showing 1.71× greater anisotropy divergence than upright ecotypes.

2. **Geodesic paths show smooth transitions** in integration, suggesting evolutionary accessibility of intermediate G-matrix structures.

3. **Tensor and affine frameworks are complementary**: Tensors identify *what* diverged (architectural vs. leaf traits); affine methods quantify *how* divergence partitioned into scale, anisotropy, and orientation.

4. **High integration facilitated rapid morphological evolution** by aligning genetic variance with the direction of selection on plant architecture.

These results demonstrate that G-matrix evolution during adaptive radiation involves both changes in trait correlation structure (well-captured by tensor eigenvectors) and changes in integration strength (uniquely quantified by affine anisotropy). The framework provides a geometrically principled approach to measuring and interpreting G-matrix divergence that complements existing multivariate quantitative genetic methods.

---

## References

Aguirre, J. D., Hine, E., McGuigan, K., & Blows, M. W. (2014). Comparing G: multivariate analysis of genetic variation in multiple populations. *Heredity*, 112, 21–29.

Hansen, T. F., & Houle, D. (2008). Measuring and comparing evolvability and constraint in multivariate characters. *Journal of Evolutionary Biology*, 21, 1201–1219.

Roda, F., Ambrose, L., Walter, G. M., Liu, H. L., Schaul, A., Lowe, A., et al. (2013). Genomic evidence for the parallel evolution of coastal forms in the Senecio lautus complex. *Molecular Ecology*, 22, 2941–2952.

Walter, G. M., Aguirre, J. D., Blows, M. W., & Ortiz-Barrientos, D. (2018). Evolution of genetic variance during adaptive radiation. *American Naturalist*, 191, E108–E128.
