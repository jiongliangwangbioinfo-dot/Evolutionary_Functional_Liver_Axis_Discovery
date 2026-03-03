# Functional Potential Field Model Spec (Anchor-to-Axis Version)

## 1. Inputs and symbols

- `x_i ∈ R^2 / R^3`: spot/cell spatial coordinate
- `r_i`: RNA feature vector
- `m_i`: metabolite feature vector
- `h_i` (optional): morphology feature vector
- `z_i`: shared multi-modal embedding
- `u_i`: latent functional potential at `x_i`

## 2. Spatial graph

Build neighborhood graph `G=(V,E)` via kNN/radius:

- `w_ij = exp(-||x_i-x_j||^2 / sigma_x^2)`

## 3. Vascular-like anchor discovery (marker-free)

Identify anchor candidates `A={a_k}` using fused evidence:

- topology score `T_i` (centrality / channel-likeness / boundary relation)
- multi-omics transition score `M_i` (local gradient breakpoint in `z_i`)
- optional morphology score `H_i`

Combined score:

`S_i = αT_i + βM_i + γH_i`

Select local maxima of `S_i` under non-maximum suppression; keep confidence `p_i`.

## 4. Local axis construction between neighboring anchors

For neighboring anchors `(a_p, a_q)`:

- define path `P_pq` on graph by geodesic/diffusion shortest path
- define axis pseudo-coordinate `s_pq(v) ∈ [0,1]` for nodes on/near `P_pq`

Collect local axes:

- `Ξ = {ξ_t}`

## 5. Per-axis molecular profile

For each axis `ξ_t`:

- smooth expression curves along `s_t`:
  - `g_{k,t}(s)` for gene `k`
  - `b_{l,t}(s)` for metabolite `l`
- derive axis feature vector `f_t`:
  - monotonicity score
  - effect size (end-to-end delta)
  - nonlinearity
  - RNA-metabolite concordance

## 6. Axis clustering and polarity alignment

Cluster `{f_t}` to obtain axis archetypes:

- strong monotonic
- multi-peak
- weak/noisy

Polarity alignment rule:

- choose direction maximizing cross-modal consistency
- tie-break by agreement with global potential gradient `∇u`

Output aligned axes `{\hat{ξ}_t}`.

## 7. Global potential field reconstruction

Joint objective:

`L = L_align + λ_smooth L_smooth + λ_axis L_axis + λ_reg L_reg`

Where:

- `L_align`: RNA/metabolite/(morphology) representation agreement
- `L_smooth = Σ_(i,j∈E) w_ij (u_i-u_j)^2`
- `L_axis`: consistency between local aligned axes and potential ordering
- `L_reg`: anti-collapse / sparsity / stability regularization

Normalized comparable pseudo-space:

`u*_i = (u_i - min(u)) / (max(u)-min(u))`

## 8. Pattern classification

Classify sample into:

- `Global axis`
- `Multi-local axes`
- `Weak organization`
- `Absent`

Using statistics:

- `S_global`: dominant gradient directionality
- `N_local`: effective number of stable local axis units
- `C_modal`: RNA-metabolite consistency
- permutation p-value against spatially shuffled baseline

## 9. PV-like / CV-like endpoint inference

Infer endpoint roles functionally (not anatomy labels):

- PV-like: low-potential anchor neighborhoods
- CV-like: high-potential anchor neighborhoods

with confidence intervals from bootstrap/permutation.
