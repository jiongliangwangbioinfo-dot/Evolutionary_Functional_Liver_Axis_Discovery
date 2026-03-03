#!/usr/bin/env python3
"""Generate simulated spatial multi-omics datasets for testing prototype code."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import mean

from prototype import simulate


def write_dataset(outdir: Path, prefix: str, n_side: int, seed: int) -> dict:
    x, rna, met, u_true = simulate(n_side=n_side, seed=seed)
    n_spots = len(x)
    n_genes = len(rna[0]) if rna else 0
    n_mets = len(met[0]) if met else 0

    spots_path = outdir / f"{prefix}_spots.csv"
    rna_path = outdir / f"{prefix}_rna.csv"
    met_path = outdir / f"{prefix}_met.csv"

    with spots_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["spot_id", "x", "y", "u_true"])
        for i, (coord, u) in enumerate(zip(x, u_true)):
            w.writerow([i, coord[0], coord[1], u])

    with rna_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["spot_id"] + [f"gene_{i}" for i in range(n_genes)])
        for i, row in enumerate(rna):
            w.writerow([i] + row)

    with met_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["spot_id"] + [f"met_{i}" for i in range(n_mets)])
        for i, row in enumerate(met):
            w.writerow([i] + row)

    summary = {
        "dataset": prefix,
        "seed": seed,
        "n_side": n_side,
        "n_spots": n_spots,
        "n_genes": n_genes,
        "n_metabolites": n_mets,
        "u_true_mean": mean(u_true),
        "u_true_min": min(u_true),
        "u_true_max": max(u_true),
        "files": {
            "spots": spots_path.name,
            "rna": rna_path.name,
            "met": met_path.name,
        },
    }
    return summary


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="data/simulated")
    p.add_argument("--n-side", type=int, default=24)
    p.add_argument("--seeds", nargs="+", type=int, default=[1, 7, 42])
    p.add_argument("--prefix", default="sim")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_summaries = []
    for idx, seed in enumerate(args.seeds):
        name = f"{args.prefix}_{idx}_seed{seed}"
        summary = write_dataset(outdir, name, args.n_side, seed)
        all_summaries.append(summary)

    summary_path = outdir / "manifest.json"
    summary_path.write_text(json.dumps(all_summaries, indent=2), encoding="utf-8")

    print(f"Generated {len(all_summaries)} dataset(s) in {outdir}")
    print(f"Manifest: {summary_path}")


if __name__ == "__main__":
    main()
