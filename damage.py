#!/usr/bin/env python3
"""Left damage-spreading speed of Rule 30: flip one cell in a random row, track the leftmost
cell where the two evolutions differ. Compared against the frozen-wedge frontier speed."""
import json
import sys

import numpy as np

T, W, SEEDS = int(sys.argv[3]) if len(sys.argv) > 3 else 20000, 0, int(sys.argv[2]) if len(sys.argv) > 2 else 16
W = 2 * T + 64
rng = np.random.default_rng(42)
speeds = []
for k in range(SEEDS):
    a = rng.integers(0, 2, W, dtype=np.uint8)
    b = a.copy()
    c = W - T - 32          # flip near the right so the left front has room
    b[c] ^= 1
    left = []
    for t in range(T):
        a[1:-1] = a[:-2] ^ (a[1:-1] | a[2:])
        b[1:-1] = b[:-2] ^ (b[1:-1] | b[2:])
        d = np.flatnonzero(a != b)
        left.append(c - d[0] if len(d) else 0)
    left = np.array(left, dtype=float)
    t = np.arange(1, T + 1)
    sl = np.polyfit(t[T // 4:], left[T // 4:], 1)[0]
    speeds.append(sl)
    print(f"seed {k}: left front speed {sl:.5f}", flush=True)
res = {"T": T, "seeds": SEEDS, "speeds": speeds, "mean": float(np.mean(speeds)),
       "sem": float(np.std(speeds, ddof=1) / np.sqrt(len(speeds)))}
json.dump(res, open(sys.argv[1], "w"), indent=1)
print(res["mean"], res["sem"])
