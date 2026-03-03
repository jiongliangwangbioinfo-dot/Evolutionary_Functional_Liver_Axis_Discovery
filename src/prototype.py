#!/usr/bin/env python3
"""Minimal runnable scaffold for anchor-to-axis PV-CV functional field recovery.

No external dependencies (stdlib only).
"""

from __future__ import annotations

import argparse
import heapq
import math
import random
from statistics import mean


def dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)


def corr(x, y):
    n = len(x)
    if n == 0:
        return 0.0
    mx, my = sum(x) / n, sum(y) / n
    vx = sum((v - mx) ** 2 for v in x)
    vy = sum((v - my) ** 2 for v in y)
    if vx <= 1e-12 or vy <= 1e-12:
        return 0.0
    c = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    return c / math.sqrt(vx * vy)


def normalize(v):
    mn, mx = min(v), max(v)
    if mx - mn <= 1e-12:
        return [0.0 for _ in v]
    return [(x - mn) / (mx - mn) for x in v]


def simulate(n_side=30, seed=0):
    random.seed(seed)
    x = []
    for i in range(n_side):
        for j in range(n_side):
            x.append((i / (n_side - 1), j / (n_side - 1)))

    c1, c2 = (0.2, 0.3), (0.75, 0.7)
    u_true = []
    for p in x:
        v = math.exp(-(dist(p, c2) ** 2) / 0.04) - 0.9 * math.exp(-(dist(p, c1) ** 2) / 0.03) + 0.15 * p[0]
        u_true.append(v)
    u_true = normalize(u_true)

    n_genes, n_mets = 24, 12
    wg = [random.gauss(1.0, 0.3) for _ in range(n_genes)]
    wm = [random.gauss(0.9, 0.25) for _ in range(n_mets)]

    rna, met = [], []
    for u in u_true:
        rg = []
        for k in range(n_genes):
            rg.append(u * wg[k] + 0.2 * (u**2) * random.gauss(1, 0.2) + random.gauss(0, 0.08))
        mt = []
        for k in range(n_mets):
            mt.append(u * wm[k] + 0.15 * math.sin(math.pi * u) * random.gauss(1, 0.2) + random.gauss(0, 0.07))
        rna.append(rg)
        met.append(mt)
    return x, rna, met, u_true


def build_knn_graph(x, k=8):
    n = len(x)
    neighbors = [[] for _ in range(n)]
    weights = [[] for _ in range(n)]
    for i in range(n):
        ds = sorted([(dist(x[i], x[j]), j) for j in range(n) if j != i], key=lambda t: t[0])[:k]
        med = ds[len(ds) // 2][0] if ds else 1.0
        for d, j in ds:
            neighbors[i].append(j)
            weights[i].append(math.exp(-(d * d) / (med * med + 1e-8)))
    return neighbors, weights


def transition_score(z, neighbors):
    n = len(z)
    sc = [0.0] * n
    for i in range(n):
        vals = []
        for j in neighbors[i]:
            dx = z[j][0] - z[i][0]
            dy = z[j][1] - z[i][1]
            vals.append(math.sqrt(dx * dx + dy * dy))
        sc[i] = sum(vals) / len(vals) if vals else 0.0
    return normalize(sc)


def detect_anchors(x, z, neighbors, n_anchors=24, min_dist=0.06):
    s = transition_score(z, neighbors)
    order = sorted(range(len(s)), key=lambda i: -s[i])
    out = []
    for idx in order:
        if len(out) >= n_anchors:
            break
        if all(dist(x[idx], x[j]) >= min_dist for j in out):
            out.append(idx)
    return out


def shortest_path(x, neighbors, src, dst):
    n = len(x)
    d = [float("inf")] * n
    prev = [-1] * n
    d[src] = 0.0
    pq = [(0.0, src)]
    while pq:
        cd, u = heapq.heappop(pq)
        if cd > d[u]:
            continue
        if u == dst:
            break
        for v in neighbors[u]:
            nd = cd + dist(x[u], x[v])
            if nd < d[v]:
                d[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))
    if d[dst] == float("inf"):
        return []
    path = [dst]
    cur = dst
    while cur != src:
        cur = prev[cur]
        if cur == -1:
            return []
        path.append(cur)
    path.reverse()
    return path


def build_axes(x, neighbors, anchors, max_pairs=40):
    pairs = {}
    for i, a in enumerate(anchors):
        dlist = sorted([(dist(x[a], x[b]), b) for j, b in enumerate(anchors) if j != i], key=lambda t: t[0])[:3]
        for d, b in dlist:
            key = tuple(sorted((a, b)))
            pairs[key] = min(d, pairs.get(key, float("inf")))
    axes = []
    for (a, b), _ in sorted(pairs.items(), key=lambda t: t[1])[:max_pairs]:
        p = shortest_path(x, neighbors, a, b)
        if len(p) >= 4:
            axes.append(p)
    return axes


def mean_profile(mat, axis, bins=10):
    prof = []
    n = len(axis)
    for b in range(bins):
        lo = b / bins
        hi = (b + 1) / bins
        idx = []
        for t, node in enumerate(axis):
            s = t / (n - 1) if n > 1 else 0
            if lo <= s <= hi + 1e-12:
                idx.append(node)
        if not idx:
            prof.append(0.0)
            continue
        vals = [mean(mat[i]) for i in idx]
        prof.append(mean(vals))
    return prof


def align_axis(rna, met, axis, proxy_u):
    p1, q1 = mean_profile(rna, axis), mean_profile(met, axis)
    s1 = corr(p1, q1)
    rev = list(reversed(axis))
    p2, q2 = mean_profile(rna, rev), mean_profile(met, rev)
    s2 = corr(p2, q2)

    # polarity: keep endpoint with larger proxy potential as CV-like endpoint
    a, b = axis[0], axis[-1]
    prefer_rev = proxy_u[a] > proxy_u[b]

    if prefer_rev:
        return rev, max(s1, s2)
    return axis, max(s1, s2)


def regress_u(n, neighbors, weights, constraints, iters=250):
    u = [0.5] * n
    fixed = {}
    for i, v in constraints:
        fixed[i] = v
        u[i] = v

    for _ in range(iters):
        nu = u[:]
        for i in range(n):
            if i in fixed:
                nu[i] = fixed[i]
                continue
            den = sum(weights[i])
            if den <= 1e-12:
                continue
            num = 0.0
            for w, j in zip(weights[i], neighbors[i]):
                num += w * u[j]
            nu[i] = num / den
        u = nu
    return normalize(u)


def classify(u, x, n_axes):
    gx = abs(corr(u, [p[0] for p in x]))
    gy = abs(corr(u, [p[1] for p in x]))
    gs = max(gx, gy)
    if gs > 0.75 and n_axes < 20:
        return "Global axis"
    if n_axes >= 10:
        return "Multi-local axes"
    if (max(u) - min(u)) > 0.1:
        return "Weak organization"
    return "Absent"


def run(n_side=30, seed=0, k=8, n_anchors=24):
    x, rna, met, u_true = simulate(n_side=n_side, seed=seed)
    neighbors, weights = build_knn_graph(x, k=k)
    z = [(mean(rna[i]), mean(met[i])) for i in range(len(x))]
    anchors = detect_anchors(x, z, neighbors, n_anchors=n_anchors)
    axes = build_axes(x, neighbors, anchors)

    proxy_u = normalize([0.5 * (mean(rna[i]) + mean(met[i])) for i in range(len(x))])

    aligned, scores = [], []
    for ax in axes:
        aa, sc = align_axis(rna, met, ax, proxy_u)
        aligned.append(aa)
        scores.append(sc)

    constraints = []
    for ax in aligned:
        constraints.append((ax[0], 0.0))
        constraints.append((ax[-1], 1.0))
    u_hat = regress_u(len(x), neighbors, weights, constraints)

    print("=== Minimal Anchor-to-Axis Prototype (Simulated Data) ===")
    print(f"spots: {len(x)}")
    print(f"anchors detected: {len(anchors)}")
    print(f"local axes built: {len(axes)}")
    print(f"mean axis alignment score: {mean(scores) if scores else 0.0:.4f}")
    print(f"corr(u_true, u_hat): {corr(u_true, u_hat):.4f}")
    print(f"classification (heuristic): {classify(u_hat, x, len(aligned))}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n-side", type=int, default=30)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--k", type=int, default=8)
    p.add_argument("--n-anchors", type=int, default=24)
    a = p.parse_args()
    run(n_side=a.n_side, seed=a.seed, k=a.k, n_anchors=a.n_anchors)


if __name__ == "__main__":
    main()
