diff --git a/AA.md b/AA.md
index 7245377826f843db4b6349cfb0ce6d0f7a90bd74..f8507435d682e1b51201b75ada24a294c5e5d136 100644
--- a/AA.md
+++ b/AA.md
@@ -1,273 +1,159 @@
-# Cross-species Functional Liver Axis Discovery Framework
+# Cross-species Functional Perfusion Field Discovery (PV–CV Conserved, Marker-free)
 
-## 1. Objective
+## 1) Biological premise and research target
 
-Given spatial transcriptomics, spatial metabolomics, coordinates, and optional histology features for each species, infer whether a **dominant functional spatial axis** exists **without predefined PV/CV landmarks**. If present, quantify its direction, strength, and cross-modal support, then align species into a shared pseudo-space `u \in [0,1]`.
+We study vertebrate liver spatial multi-omics with:
+- spatial transcriptomics,
+- spatial metabolomics,
+- spatial coordinates,
+- optional morphology.
 
----
-
-## 2. Inputs and Notation
-
-For species `s`:
+Assumption:
+- PV-like inflow and CV-like outflow circulation exists across vertebrates.
 
-- Spots/cells: `i = 1...N_s` with coordinates `x_i \in R^2` (or `R^3`)
-- RNA matrix: `G_s \in R^{N_s \times P_s}`
-- Metabolite matrix: `M_s \in R^{N_s \times Q_s}`
-- Optional morphology features: `H_s \in R^{N_s \times R_s}` (texture/CNN embeddings from H&E)
-- Optional phylogenetic metadata: evolutionary distance/time for downstream trend testing
+No assumption:
+- no fixed mammalian lobule template,
+- no fixed PV/CV geometry,
+- no guaranteed single global zonation axis.
 
-No assumption is made about vessel labels, triads, or marker conservation.
+Core question:
+- whether conserved perfusion is functionally spatialized into detectable organization.
 
 ---
 
-## 3. High-level Pipeline
-
-1. **Within-species preprocessing (modality-specific denoising + normalization)**
-2. **Spatial graph construction**
-3. **Unsupervised multi-modal latent embedding** with spatial smoothness
-4. **Dominant gradient/axis extraction** from latent field
-5. **Axis existence and strength testing** against spatial nulls
-6. **Cross-modal consistency scoring**
-7. **Axis orientation and pseudo-space assignment `u \in [0,1]`**
-8. **Cross-species alignment and evolutionary analysis**
-
----
-
-## 4. Step-by-step Algorithm
-
-## 4.1 Preprocess each modality
-
-For each species `s`:
-
-- RNA:
-  - Library-size normalization, variance stabilization (e.g., SCTransform/log1p)
-  - Remove low-quality spots/features
-  - Select highly variable genes
-- Metabolites:
-  - TIC/quantile normalization, log transform, batch correction if needed
-- Morphology (optional):
-  - Extract patch-level descriptors (e.g., self-supervised CNN features)
-  - Standardize features
-
-Result: normalized matrices `G'_s, M'_s, H'_s`.
-
-## 4.2 Build spatial neighborhood graph
-
-Construct `W_s` over coordinates via kNN or radius graph, optionally anisotropic if tissue geometry demands.
-
-- Graph Laplacian `L_s = D_s - W_s`
-- Preserve local continuity while allowing global gradient detection
-
-## 4.3 Learn a shared latent spatial representation
-
-Learn latent vector `z_i \in R^d` for each spot using a multi-view objective:
-
-`Loss = L_reconRNA + L_reconMet + (optional) L_reconMorph + lambda * Tr(Z^T L_s Z) + beta * L_align`
-
-Where:
-
-- reconstruction losses are modality-specific (NB/Gaussian/Poisson as appropriate)
-- `Tr(Z^T L_s Z)` enforces spatial smoothness
-- `L_align` encourages shared structure across modalities (CCA/contrastive co-embedding)
-
-Implementation options (choose one):
-
-- Multi-omic variational graph autoencoder
-- MOFA+ with spatial regularization
-- Graph-regularized joint NMF/CCA
-
-Output per species: latent matrix `Z_s \in R^{N_s \times d}`.
-
-## 4.4 Extract candidate dominant axis
-
-Goal: find scalar field `f_i` (pseudo-coordinate before scaling) that captures strongest smooth spatial variation.
-
-Two equivalent approaches:
-
-1. **Principal smooth component**:
-   - Solve generalized eigenproblem maximizing variance under smoothness prior.
-2. **Diffusion geometry**:
-   - Compute diffusion map on `Z_s` + spatial kernel
-   - Use first non-trivial diffusion component `psi_1` as candidate axis
-
-Set `f = psi_1` (or top smooth component).
-
-Convert to unit interval:
-
-`u_i = (f_i - min(f)) / (max(f)-min(f)) \in [0,1]`
-
-`u` is the species-specific pseudo-space.
-
-## 4.5 Decide whether a dominant axis exists
-
-Define **axis strength** from three complementary criteria:
-
-1. **Explained smooth variance**
-   - `S_var = lambda1 / sum_k lambda_k` from smooth-spectrum decomposition
-2. **Anisotropy / dominance gap**
-   - `S_gap = (lambda1 - lambda2) / lambda1`
-3. **Spatial autocorrelation of `u`**
-   - Moran's I or Geary's C (`S_auto`)
-
-Aggregate strength:
-
-`S_axis = w1*S_var + w2*S_gap + w3*S_auto` (weights via CV or fixed equal weights).
+## 2) Computational objective
 
-### Statistical significance
+Infer, without vessel markers or anatomical priors:
 
-Generate null distribution via spatially constrained permutations:
+1. Whether perfusion-related functional spatial organization exists.
+2. Which organization type is present:
+   - **Global axis**: one dominant tissue-scale gradient.
+   - **Multi-local axes**: repeated local PV–CV-like units.
+   - **Weak organization**: limited gradient structure.
+   - **Absent**: no robust organization.
+3. A species-comparable pseudo-space `u*(x) in [0,1]`.
 
-- Permute feature vectors while preserving local neighborhood statistics, or
-- Randomly rotate/warp coordinates within tissue mask while preserving density
-
-Compute `p_axis = P_null(S_axis^null >= S_axis^obs)`.
-
-Decision:
-
-- If `p_axis < alpha` and `S_axis > tau` -> axis exists
-- Else -> no dominant functional axis detected
-
-This explicitly handles species lacking zonation.
-
-## 4.6 Recover axis direction in physical space
-
-Estimate gradient field in tissue coordinates:
-
-`g(x_i) = nabla_x u_i`
-
-(using local linear regression over neighbors).
-
-Outputs:
-
-- global direction vector (mean normalized gradient)
-- local vector field map (to inspect flow-like organization)
-
-Because PV/CV identity is unknown, orientation sign is arbitrary at this stage.
-
-## 4.7 Cross-modal consistency score
-
-Independently estimate pseudo-space per modality:
-
-- `u^RNA`, `u^MET`, `u^MORPH` (if available)
-
-Then compute:
-
-- Rank correlation: `rho_RM = Spearman(u^RNA, u^MET)`
-- Distance correlation / mutual information between modal coordinates
-- Procrustes similarity between 1D embeddings
+---
 
-Composite consistency:
+## 3) Modeling idea: continuous functional potential field
 
-`C_multi = a1*rho_RM + a2*rho_RH + a3*rho_MH` (using available pairs)
+Instead of forcing one linear axis, model tissue as a scalar field `u(x)`:
 
-Also test against null permutations for significance.
+- jointly supported by RNA + metabolite (+ morphology),
+- spatially smooth but allowed to be multi-basin / irregular,
+- interpretable via local gradient flow `∇u(x)`.
 
-High `C_multi` indicates conserved cross-modal gradient rather than modality-specific artifact.
+Interpretation:
+- CV-like regions: local high potential,
+- PV-like regions: local low potential,
+- local functional axis: gradient trajectories from low to high potential.
 
 ---
 
-## 5. Cross-species alignment to shared pseudo-space
+## 4) Algorithm pipeline
 
-Each species has `u_s \in [0,1]` but orientation/scale can differ.
+### 4.1 Modality preprocessing
+For each species `s`:
+- normalize and denoise RNA / metabolite / morphology,
+- keep informative features,
+- optional batch correction.
 
-## 5.1 Orientation harmonization (sign disambiguation)
+### 4.2 Spatial graph
+Build `W_s` (kNN/radius graph) from coordinates and Laplacian `L_s`.
 
-Since `u` can be flipped, orient by maximizing consistency with:
+### 4.3 Multi-modal latent fusion with spatial regularization
+Learn spot embeddings `Z_s` from all available modalities:
 
-- cross-modal agreement (choose sign maximizing `C_multi`), and/or
-- monotonicity of conserved functional programs (e.g., oxidation-to-glycolysis module scores discovered unsupervised)
+`Loss = L_RNA + L_MET + (optional) L_MORPH + λ Tr(Z_s^T L_s Z_s) + β L_cross_modal`
 
-Apply `u <- 1-u` when needed.
+Output: shared latent representation preserving cross-modal and spatial consistency.
 
-## 5.2 Nonlinear alignment across species
+### 4.4 Functional field estimation
+Estimate smooth scalar field `u(x)` from `Z_s` on the spatial graph:
+- diffusion / Laplacian smooth components,
+- graph-regularized regression / variational field models.
 
-To compare species with different gradient shapes, align distributions via:
+Normalize to `u*(x) in [0,1]` for comparability.
 
-- monotone warping `phi_s(u_s)` learned by optimal transport or quantile mapping
-- anchor features: conserved gene/metabolite modules that vary monotonically with `u`
+### 4.5 Structure typing (global vs local vs weak vs absent)
+Compute:
+- global anisotropy (dominance of first smooth component),
+- local extrema density and basin topology (Morse/Watershed on `u`),
+- gradient coherence (alignment statistics of `∇u`),
+- spatial autocorrelation.
 
-Final aligned coordinate:
+Classify:
+- **Global axis**: high dominance + low extrema density.
+- **Multi-local axes**: multiple significant extrema/basins + coherent local gradients.
+- **Weak organization**: weak but non-random structure.
+- **Absent**: not different from spatial null.
 
-`u*_s = phi_s(u_s) in [0,1]`
+### 4.6 Local functional unit recovery
+From `u(x)`:
+- identify local minima / maxima (PV-like / CV-like candidates),
+- segment basins (lobule-like functional units),
+- estimate per-unit principal flow direction via `∇u`.
 
-Now species are comparable in a shared pseudo-space.
+### 4.7 Cross-modal consistency validation
+Infer modality-specific fields `u_RNA`, `u_MET`, (`u_MORPH`) and evaluate:
+- rank correlation / mutual information,
+- basin overlap and gradient-direction agreement,
+- permutation-based significance.
 
----
+A structure is considered biologically supported only if multi-modal consistency is significant.
 
-## 6. Evolutionary inference outputs
+### 4.8 Cross-species alignment
+Because geometry may differ across species:
+- do not align by anatomical templates,
+- align by monotone warping of `u*` and conserved functional programs.
 
-For each species report:
+Result: comparable coordinate `u*_s(x) in [0,1]` for evolutionary analysis.
 
-1. **Axis existence:** boolean + p-value
-2. **Axis strength:** `S_axis` (with components)
-3. **Direction:** global vector + local vector field
-4. **Cross-modal consistency:** `C_multi` + p-value
-5. **Pseudo-space coordinate:** `u_i` (and aligned `u*_i`)
+---
 
-Across species test trends (phylogenetically aware regression):
+## 5) Evolution-oriented outputs
 
-- `S_axis ~ evolutionary_time`
-- `C_multi ~ evolutionary_time`
+Per species:
+1. `organization_type` ∈ {global, multi_local, weak, absent}
+2. global-axis strength score
+3. number of local axes / extrema / basins
+4. gradient coherence metrics
+5. cross-modal consistency score and p-value
+6. spot-level `u*(x)` and local vector field
 
-to quantify emergence/strengthening of zonation through evolution.
+Across species:
+- test whether functional spatialization strengthens through evolution,
+- quantify when conserved PV–CV circulation becomes functionally zonated.
 
 ---
 
-## 7. Practical pseudocode
+## 6) Minimal pseudocode
 
 ```text
-for species s in Species:
-    G', M', H' = preprocess_modalities(G, M, H)
-    W, L = build_spatial_graph(coords)
-
-    Z = fit_multimodal_spatial_latent(G', M', H', L)
-
-    f = first_smooth_diffusion_component(Z, W)
-    u = minmax(f)
-
-    S_var, S_gap, S_auto = compute_axis_metrics(u, Z, W)
-    S_axis = aggregate(S_var, S_gap, S_auto)
-    p_axis = spatial_null_test(S_axis, G', M', H', coords)
-    axis_exists = (p_axis < alpha) and (S_axis > tau)
-
-    grad_field, global_dir = estimate_gradient_field(u, coords, W)
-
-    u_RNA, u_MET, u_MORPH = modality_specific_axes(...)
-    C_multi, p_multi = cross_modal_consistency(u_RNA, u_MET, u_MORPH)
-
-    u_oriented = resolve_sign(u, criterion=C_multi+conserved_modules)
-    save species outputs
-
-u_star_all = cross_species_monotone_alignment({u_oriented_s}, anchors=conserved_modules)
-run_evolutionary_trend_models(S_axis, C_multi, phylogeny)
+for species s:
+    preprocess modalities
+    build spatial graph W_s, L_s
+    learn multimodal latent Z_s with spatial smoothness
+    estimate functional field u_s(x)
+    normalize to u*_s in [0,1]
+
+    compute global anisotropy + local topology + gradient coherence
+    run spatial null tests
+    classify organization_type
+
+    recover local extrema / basins / local flow axes
+    compute cross-modal consistency and significance
+
+align {u*_s} across species by monotone functional alignment
+run evolutionary trend models on organization metrics
 ```
 
 ---
 
-## 8. Why this satisfies constraints
-
-- **No PV/CV priors:** fully unsupervised gradient extraction
-- **No bile-duct/triad assumptions:** morphology optional, not required
-- **Architecture-agnostic:** graph + diffusion handles disordered tissues
-- **Can detect absence:** explicit null-tested axis existence decision
-- **Multi-modal robustness:** consistency score prevents single-modality over-interpretation
-- **Cross-species comparable:** aligned pseudo-space `[0,1]` supports evolutionary analysis
-
----
-
-## 9. Suggested diagnostics and failure modes
-
-- If `S_axis` high but `C_multi` low: possible modality artifact or spatial resolution mismatch
-- If all species non-significant: may indicate weak zonation or underpowered sampling
-- If `u` correlates strongly with section edge artifacts: include tissue-mask covariates and edge-aware nulls
-- If multiple comparable eigencomponents: tissue may contain multiple competing gradients (branching organization)
-
----
-
-## 10. Minimal deliverable table per species
-
-| Species | Axis exists | S_axis | p_axis | C_multi | p_multi | Global direction | Notes |
-|---|---:|---:|---:|---:|---:|---|---|
-| S1 | 0/1 | value | value | value | value | (dx,dy) | comments |
+## 7) Why this matches the project constraints
 
-This table plus spot-level `u`/`u*` enables all downstream biological interpretation tasks requested.
+- assumes only conserved perfusion logic, not conserved geometry,
+- explicitly supports multi-axis and irregular structures,
+- permits true “no spatial functionalization” outcome,
+- requires cross-modal support to avoid single-modality artifacts,
+- yields cross-species comparable pseudo-space for evolutionary inference.
