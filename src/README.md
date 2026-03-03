# Minimal runnable prototype (`src/`)

Run on simulated data:

```bash
python3 src/prototype.py --n-side 30 --seed 0 --k 8 --n-anchors 24
```

Pipeline implemented:

1. Anchor detection (marker-free, transition score on multi-omics embedding)
2. Local axis construction between neighboring anchors (graph geodesic path)
3. Axis polarity alignment (maximize RNA-metabolite profile agreement)
4. Global functional potential regression `u(x)` with axis endpoint constraints

Outputs:

- number of anchors
- number of local axes
- mean alignment score
- correlation between simulated `u_true` and estimated `u_hat`
- heuristic pattern classification
