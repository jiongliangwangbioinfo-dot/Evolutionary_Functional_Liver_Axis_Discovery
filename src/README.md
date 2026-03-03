# Minimal runnable prototype (`src/`)

## 1) Run anchor-to-axis prototype on simulated data

```bash
python3 src/prototype.py --n-side 30 --seed 0 --k 8 --n-anchors 24
```

Pipeline implemented:

1. Anchor detection (marker-free, transition score on multi-omics embedding)
2. Local axis construction between neighboring anchors (graph geodesic path)
3. Axis polarity alignment (cross-modal profile + proxy potential)
4. Global functional potential regression `u(x)` with axis endpoint constraints

Outputs:

- number of anchors
- number of local axes
- mean alignment score
- correlation between simulated `u_true` and estimated `u_hat`
- heuristic pattern classification

## 2) Generate reusable simulated datasets for testing

```bash
python3 src/generate_sim_data.py --outdir data/simulated --n-side 24 --seeds 1 7 42
```

Generated files (per seed):

- `*_spots.csv`: `spot_id, x, y, u_true`
- `*_rna.csv`: spot-level gene matrix
- `*_met.csv`: spot-level metabolite matrix
- `manifest.json`: dataset summary and file index
