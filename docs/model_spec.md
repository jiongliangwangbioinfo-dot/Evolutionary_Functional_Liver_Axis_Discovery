# Functional Potential Field Model Spec

## Symbols

- `x_i ∈ R^2 / R^3`: spatial coordinate of spot/cell `i`
- `r_i`: RNA feature vector
- `m_i`: metabolite feature vector
- `h_i` (optional): morphology feature vector
- `u_i`: latent functional potential at `x_i`

## Graph construction

Build spatial neighborhood graph `G=(V,E)` with kNN or radius edges.

- Edge weight `w_ij = exp(-||x_i-x_j||^2 / sigma_x^2)`.

## Objective

Learn shared representation and scalar potential jointly:

`L = L_align + λ_smooth L_smooth + λ_cons L_cons + λ_reg L_reg`

Where:

- `L_align`: multi-modal representation agreement (RNA/metabolite/(morphology))
- `L_smooth = Σ_(i,j∈E) w_ij (u_i-u_j)^2`
- `L_cons`: modality-specific potential consistency term
- `L_reg`: regularization (e.g., variance floor for non-collapsed field)

## Derived quantities

- Local gradient proxy at node `i`:
  - `g_i = Σ_(j∈N(i)) w_ij (u_j-u_i) * (x_j-x_i) / ||x_j-x_i||`
- Local extrema:
  - maxima if `u_i >= u_j, ∀j∈N(i)` (CV-like)
  - minima if `u_i <= u_j, ∀j∈N(i)` (PV-like)
- Basin segmentation:
  - assign each node by gradient ascent destination.

## Pattern classification heuristics

- **Global axis**: high first principal direction dominance of `g_i`
- **Multi-local axes**: multiple stable extrema and basins with moderate within-basin coherence
- **Weak organization**: low gradient magnitude but non-random smoothness
- **Absent**: no significant smooth structure vs permutation baseline

## Cross-species comparable outputs

- normalized pseudo-space: `u*_i = (u_i - min(u)) / (max(u)-min(u))`
- summary metrics:
  - `S_global`: global axis strength
  - `N_local`: effective number of local axis units
  - `C_modal`: cross-modal consistency score
