#!/usr/bin/env python3
"""X-File 6: is Rule 30's 2-D picture encoding a hidden third dimension?

Four signatures a hidden coordinate would leave, each compared against Rule 30 run from
fair random bits (same rule, no special initial condition) and split by region:
  T1 patch information  -- entropy of k x k patches vs the light-cone bound 3k-2 bits
  T2 row correlations   -- corr(x_s, x_s+D) within rows, D = 1..2048
  T3 scale index        -- interior statistics conditioned on v2(n+1), the 2-adic "scale"
  T4 renormalisation    -- predictability of majority-coarse-grained cells from the coarse row above
Writes thirddim.json to the output directory.
"""
import json
import sys

import numpy as np

OUT = sys.argv[1]
N = 16384
rng = np.random.default_rng(30)


def evolve_single(N):
    W = 2 * N + 1
    g = np.zeros((N, W), dtype=np.uint8)
    row = np.zeros(W + 2, dtype=np.uint8)
    row[N + 1] = 1
    for n in range(N):
        g[n] = row[1:-1]
        row[1:-1] = row[:-2] ^ (row[1:-1] | row[2:])
    return g


def evolve_ring(N, W):
    g = np.zeros((N, W), dtype=np.uint8)
    row = rng.integers(0, 2, W, dtype=np.uint8)
    for n in range(N):
        g[n] = row
        row = np.roll(row, 1) ^ (row | np.roll(row, -1))
    return g


def H(codes):
    _, c = np.unique(codes, return_counts=True)
    p = c / c.sum()
    return float(-(p * np.log2(p)).sum())


def region_cols(n, region):
    # column indices into the single-cell grid (centre at N) for row n
    if region == "interior":   # right of the wedge, away from the right edge
        lo, hi = -int(0.2 * n), int(0.9 * n)
    else:                      # frozen wedge: left edge + 100 .. s = -0.3 n
        lo, hi = -n + 100, -int(0.3 * n)
    return N + lo, N + hi


def patches(g, rows, colrange, k):
    codes = []
    for n in rows:
        lo, hi = colrange(n)
        block = g[n:n + k, lo:hi]
        if block.shape[0] < k:
            continue
        c = np.zeros(hi - lo - k + 1, dtype=np.int64)
        for i in range(k):
            for j in range(k):
                c = (c << 1) | block[i, j:j + hi - lo - k + 1]
        codes.append(c)
    return np.concatenate(codes)


def corr_profile(g, rows, colrange, maxd=2048):
    num = np.zeros(maxd + 1)
    cnt = np.zeros(maxd + 1)
    for n in rows:
        lo, hi = colrange(n)
        x = g[n, lo:hi].astype(np.float64) - 0.5
        f = np.fft.rfft(x, 2 * len(x))
        ac = np.fft.irfft(f * np.conj(f))[:maxd + 1]
        num += ac
        cnt += np.maximum(len(x) - np.arange(maxd + 1), 0)
    c = num / cnt / 0.25
    return c


def coarse(g, rows, colrange, k):
    """Majority of k x k blocks; returns conditional entropy H(C(t+1,j) | C(t,j-1..j+1))."""
    ctx, nxt = [], []
    for n in rows:
        lo, hi = colrange(n)
        w = (hi - lo) // k * k
        if n + 2 * k > g.shape[0]:
            continue
        a = g[n:n + k, lo:lo + w].reshape(k, w // k, k).sum(axis=(0, 2))
        b = g[n + k:n + 2 * k, lo:lo + w].reshape(k, w // k, k).sum(axis=(0, 2))
        ca = (2 * a > k * k).astype(np.int64)
        cb = (2 * b > k * k).astype(np.int64)
        ctx.append(ca[:-2] * 4 + ca[1:-1] * 2 + ca[2:])
        nxt.append(cb[1:-1])
    ctx, nxt = np.concatenate(ctx), np.concatenate(nxt)
    return H(ctx * 2 + nxt) - H(ctx)


def main():
    single = evolve_single(N)
    ring = evolve_ring(N, 2 * N)
    rows = range(N // 2, N - 16, 8)
    ring_cols = lambda n: (0, 2 * N)
    res = {"rows": [N // 2, N - 16], "N": N}

    regions = {"interior": lambda n: region_cols(n, "interior"),
               "wedge": lambda n: region_cols(n, "wedge"),
               "random_start": ring_cols}
    grids = {"interior": single, "wedge": single, "random_start": ring}

    res["T1_patch_entropy"] = {}
    for k in (2, 3, 4):
        res["T1_patch_entropy"][k] = {"bound_3k-2": 3 * k - 2, "max_k2": k * k}
        for name, cr in regions.items():
            res["T1_patch_entropy"][k][name] = H(patches(grids[name], rows, cr, k))

    res["T2_corr"] = {}
    for name, cr in regions.items():
        c = corr_profile(grids[name], range(N // 2, N - 16, 64), cr)
        res["T2_corr"][name] = {"max_abs_D1_2048": float(np.max(np.abs(c[1:]))),
                                "argmax": int(np.argmax(np.abs(c[1:])) + 1),
                                "first8": [round(float(x), 5) for x in c[1:9]],
                                "rms_D1_2048": float(np.sqrt(np.mean(c[1:] ** 2)))}

    # T3: interior density and 111-run density by v2(n+1)
    v2 = lambda x: (x & -x).bit_length() - 1
    t3 = {}
    for n in range(N // 2, N):
        lo, hi = region_cols(n, "interior")
        x = single[n, lo:hi]
        v = min(v2(n + 1), 6)
        d = t3.setdefault(v, [0, 0, 0])
        d[0] += int(x.sum())
        d[1] += len(x)
        d[2] += int((x[:-2] & x[1:-1] & x[2:]).sum())
    res["T3_scale"] = {str(v): {"rows_weighted_cells": d[1], "density": d[0] / d[1],
                                "z_density": (d[0] - d[1] / 2) / (np.sqrt(d[1]) / 2),
                                "frac111": d[2] / d[1]} for v, d in sorted(t3.items())}

    res["T4_coarse_cond_entropy"] = {}
    for k in (2, 4, 8):
        res["T4_coarse_cond_entropy"][k] = {name: coarse(grids[name], range(N // 2, N - 32, 16), cr, k)
                                            for name, cr in regions.items()}

    with open(f"{OUT}/thirddim.json", "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
