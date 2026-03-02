 (cd "$(git rev-parse --show-toplevel)" && git apply --3way <<'EOF' 
diff --git a/CROSS_SPECIES_PV_CV_AXIS_ALGORITHM.md b/CROSS_SPECIES_PV_CV_AXIS_ALGORITHM.md
new file mode 100644
index 0000000000000000000000000000000000000000..877c4af57e0fc1e6a94973db3ba5463d4cc35f85
--- /dev/null
+++ b/CROSS_SPECIES_PV_CV_AXIS_ALGORITHM.md
@@ -0,0 +1,209 @@
+# Cross-species Functional PV–CV Axis Discovery Framework
+
+## 1. Objective
+Given spatial transcriptomics, spatial metabolomics, coordinates, and optional histology-derived morphology from each species, infer whether a **dominant functional spatial axis** exists without using predefined vascular landmarks, and align species into a shared pseudo-space \([0,1]\).
+
+---
+
+## 2. Inputs and Notation
+For each species \(s\):
+
+- Coordinates: \(X_s \in \mathbb{R}^{n_s \times d}\), typically \(d=2\)
+- Gene matrix: \(G_s \in \mathbb{R}^{n_s \times p_s}\)
+- Metabolite matrix: \(M_s \in \mathbb{R}^{n_s \times q_s}\)
+- Optional morphology matrix: \(H_s \in \mathbb{R}^{n_s \times r_s}\)
+
+where \(n_s\) is number of spots/cells/pixels.
+
+---
+
+## 3. Species-wise Axis Discovery
+
+### Step 3.1: Preprocessing and graph construction
+1. Normalize each modality (e.g., SCTransform/logCPM for RNA, TIC+log for metabolites, z-score for morphology).
+2. Remove low-quality spots and low-information features.
+3. Build a spatial neighborhood graph \(\mathcal{G}_s=(V,E)\) from coordinates (kNN or radius graph).
+4. Compute graph Laplacian \(L_s\).
+
+Why graph-first: this avoids assuming any specific liver anatomy and allows irregular tissue geometry.
+
+### Step 3.2: Unsupervised multimodal latent field
+Construct modality-specific denoised embeddings:
+- \(Z_s^G = f_G(G_s)\)
+- \(Z_s^M = f_M(M_s)\)
+- \(Z_s^H = f_H(H_s)\) (if available)
+
+Use PCA/NMF/autoencoder per modality, then fuse into a shared latent representation \(Z_s\) with a graph-regularized objective:
+
+\[
+\min_{Z_s,\{W_m\}} \sum_{m\in\{G,M,H\}} \alpha_m\|X_s^mW_m - Z_s\|_F^2 + \lambda\operatorname{tr}(Z_s^TL_sZ_s)
+\]
+
+- \(X_s^m\): matrix for modality \(m\)
+- \(\alpha_m\): modality reliability weights (learned or quality-based)
+- Smoothness term enforces spatially coherent functional patterns.
+
+### Step 3.3: Candidate gradient extraction (marker-free)
+Estimate dominant smooth gradients from \(Z_s\):
+
+1. Compute diffusion operator on \(\mathcal{G}_s\) using \(Z_s\)-aware affinities.
+2. Extract first non-trivial diffusion components \(\phi_{s,1},\phi_{s,2},...\).
+3. Each \(\phi_{s,k}\) is a candidate functional axis.
+
+Choose dominant axis \(a_s\) as component maximizing both:
+- high spatial autocorrelation (Moran's I / Geary's C), and
+- high cross-modal variance explained.
+
+This detects an axis if one exists; if tissue is isotropic/noisy, no component will pass significance thresholds.
+
+### Step 3.4: Axis existence test
+Define null models by spatial permutation preserving graph degree and local intensity distribution:
+
+- Permute feature values over graph-constrained randomizations.
+- Recompute top-axis score each time.
+
+Axis existence p-value:
+\[
+p_s = \Pr(\text{null score} \ge \text{observed score})
+\]
+
+Declare dominant axis present if:
+- \(p_s < \tau_p\) (e.g., 0.01), and
+- effect size \(E_s\) above threshold (e.g., standardized axis contrast).
+
+If not, return "no dominant axis".
+
+### Step 3.5: Direction and vector field
+The scalar axis \(a_s\) is sign-ambiguous. Provide:
+
+1. **Coordinate gradient field**: \(\nabla a_s(x)\) (smoothed finite differences on graph).
+2. **Principal direction vector** \(v_s\): first principal direction of \((X_s, a_s)\)-induced gradient.
+3. Sign convention for comparability:
+   - orient so that conserved oxidative→glycolytic module score increases with pseudo-space, or
+   - orient by maximal agreement with inferred perfusion proxy if available.
+
+---
+
+## 4. Quantitative Outputs per Species
+
+### 4.1 Dominant-axis strength
+Use composite strength score:
+\[
+S_s = w_1 I_s + w_2 R^2_{\text{spatial},s} + w_3 C_s
+\]
+where:
+- \(I_s\): Moran's I of \(a_s\)
+- \(R^2_{\text{spatial},s}\): variance in \(Z_s\) explained by \(a_s\)
+- \(C_s\): contrast-to-noise between axis extremes
+
+Report confidence intervals by bootstrap over spots.
+
+### 4.2 Cross-modal consistency
+Compute modality-specific axis projections \(a_s^G,a_s^M,a_s^H\) and evaluate:
+
+- pairwise Spearman correlations
+- mutual information between binned pseudo-space and modality features
+- Procrustes alignment error between modality gradient embeddings
+
+Aggregate into consistency score:
+\[
+K_s = \beta_1\rho(G,M)+\beta_2\rho(G,H)+\beta_3\rho(M,H)-\beta_4\varepsilon_{\text{proc}}
+\]
+
+---
+
+## 5. Cross-species Pseudo-space Alignment \([0,1]\)
+
+### Step 5.1: Within-species normalization
+Convert axis to monotonic pseudo-space:
+\[
+u_s = \operatorname{rank}(a_s)/(n_s-1) \in [0,1]
+\]
+(rank normalization is robust to nonlinear scaling).
+
+### Step 5.2: Functional anchoring without anatomical landmarks
+Construct conserved functional programs (data-driven):
+1. Identify cross-species ortholog-mapped genes/metabolites.
+2. Derive shared latent factors (e.g., multi-species NMF/canonical correlation).
+3. Use factors with monotonic behavior along \(\nu_s\) as anchors.
+
+### Step 5.3: Manifold alignment across species
+Align \(\nu_s\) using monotone warping functions \(T_s:[0,1]\to[0,1]\):
+
+\[
+\min_{\{T_s\}}\sum_{s}\sum_{f\in\mathcal{F}}\|g_{s,f}(\nu)-\bar g_f(T_s(\nu))\|^2 + \gamma\|T_s''\|^2
+\]
+subject to \(T_s'\ge 0\).
+
+Output aligned coordinate:
+\[
+\tilde\nu_s = T_s(\nu_s)
+\]
+
+This supports species comparison even if architecture differs.
+
+---
+
+## 6. Evolutionary Analysis Layer
+Using \(S_s\), \(K_s\), and aligned pseudo-space \(\tilde\nu_s\):
+
+1. Test trend with phylogenetic depth:
+   - phylogenetic generalized least squares (PGLS): \(S_s \sim \text{evolutionary position}\)
+2. Quantify emergence of structured zonation:
+   - increasing \(S_s\), decreasing alignment residuals, increasing cross-modal coherence.
+3. Compare slope/shape of conserved metabolic programs along \(\tilde\nu\).
+
+Interpretation:
+- low \(S_s\), low \(K_s\): no robust axis
+- moderate \(S_s\), high \(K_s\): functional gradient present without strong anatomy
+- high \(S_s\), high \(K_s\): strongly organized zonation analogous to PV–CV function
+
+---
+
+## 7. Algorithm Pseudocode
+```text
+for species s in SpeciesSet:
+    preprocess G_s, M_s, H_s, X_s
+    build spatial graph Graph_s from X_s
+    learn shared latent Z_s via multimodal graph-regularized fusion
+    compute diffusion components phi_{s,k}
+    score each component by spatial autocorrelation + cross-modal support
+    select best component a_s
+    run graph-constrained permutation test
+    if not significant:
+        label species as "no dominant axis"
+        continue
+    orient axis sign using functional/perfusion proxy
+    compute pseudo-space nu_s in [0,1]
+    compute strength S_s and consistency K_s
+
+collect all species with valid axes
+learn monotone alignment warps T_s
+compute aligned pseudo-space tilde_nu_s
+perform evolutionary trend modeling on S_s, K_s, aligned profiles
+```
+
+---
+
+## 8. Practical Implementation Notes
+- Robustness checks:
+  - sensitivity to k in kNN graph
+  - repeated downsampling of spots
+  - leave-one-modality-out stability
+- Missing modality handling:
+  - adjust \(\alpha_m\), report reduced-confidence \(K_s\)
+- "No axis" is an explicit valid outcome, not a failure.
+
+---
+
+## 9. Deliverables (per species)
+1. `axis_exists` (boolean + p-value)
+2. `axis_scalar` (spot-level scalar field)
+3. `axis_vector_field` (local gradient vectors)
+4. `axis_strength` \(S_s\)
+5. `cross_modal_consistency` \(K_s\)
+6. `pseudo_space_raw` \(\nu_s\in[0,1]\)
+7. `pseudo_space_aligned` \(\tilde\nu_s\in[0,1]\)
+8. QC + bootstrap confidence report
+
+This framework directly targets discovery of an evolutionarily emerging, conserved **functional** zonation axis without imposing mammalian anatomy.
 
EOF
)
