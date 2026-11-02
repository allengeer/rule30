#!/usr/bin/env python3
"""X-File 1c: window-free wedge frontier vs left damage speed.
usage: frontier_analysis.py results_dir
Reads frontier_10m/frontier.i64 (L(m), distance of the leftmost period-32 break from the left
edge) and damage_*/damage_<seed>.i32 (left damage front after each step). Writes summary.json
and frontier_vs_damage.png. Speeds are leftward cells per step in absolute coordinates, so the
frontier's speed is 1 - dL/dm."""
import glob
import json
import sys

import numpy as np

R = sys.argv[1]
out = {}

L = np.fromfile(f"{R}/frontier_10m/frontier.i64", np.int64)
N = int(np.flatnonzero(L > 0)[-1]) + 1 if (L > 0).any() else 0
L = L[:N]
m = np.arange(N, dtype=float)
fr = {"rows": N, "L_over_m": {}}
for k in [10**5, 2 * 10**5, 5 * 10**5, 10**6, 2 * 10**6, 5 * 10**6, 10**7 - 1]:
    if k < N:
        fr["L_over_m"][str(k)] = float(L[k] / k)
# Global slope over the last 3/4, and slopes over independent consecutive blocks.
a = N // 4
slope, icpt = np.polyfit(m[a:], L[a:], 1)
fr["speed_fit_last3q"] = 1 - slope
fr["intercept"] = icpt
resid = L[a:] - (slope * m[a:] + icpt)
fr["resid_rms"] = float(resid.std())
blocks = []
B = 8
edges = np.linspace(a, N, B + 1).astype(int)
for i in range(B):
    s = slice(edges[i], edges[i + 1])
    blocks.append(1 - np.polyfit(m[s], L[s], 1)[0])
fr["speed_blocks"] = blocks
fr["speed_blocks_mean"] = float(np.mean(blocks))
fr["speed_blocks_sem"] = float(np.std(blocks, ddof=1) / np.sqrt(B))
# Mean front speed between far-apart rows is robust to the fluctuations: (m - L) / m at the end.
fr["speed_chord"] = float(1 - (L[N - 1] - L[a]) / (N - 1 - a))
out["frontier"] = fr

dm = {}
for sub in sorted(glob.glob(f"{R}/damage_*/")):
    files = sorted(glob.glob(f"{sub}damage_*.i32"))
    if not files:
        continue
    sp, ch = [], []
    for f in files:
        x = np.fromfile(f, np.int32).astype(float)
        T = len(x)
        t = np.arange(1, T + 1, dtype=float)
        sp.append(np.polyfit(t[T // 4:], x[T // 4:], 1)[0])
        ch.append(x[-1] / T)
    name = sub.rstrip("/").split("/")[-1]
    dm[name] = {"T": T, "seeds": len(sp), "mean": float(np.mean(sp)),
                "sem": float(np.std(sp, ddof=1) / np.sqrt(len(sp))) if len(sp) > 1 else None,
                "sd": float(np.std(sp, ddof=1)) if len(sp) > 1 else None,
                "chord_mean": float(np.mean(ch)),
                "chord_sem": float(np.std(ch, ddof=1) / np.sqrt(len(ch))) if len(ch) > 1 else None}
out["damage"] = dm
json.dump(out, open(f"{R}/summary.json", "w"), indent=1)
print(json.dumps(out, indent=1))

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit(0)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
step = max(1, N // 20000)
ax[0].plot(m[step::step], 1 - L[step::step] / m[step::step], lw=0.6, color="#1f5fa8", label="frontier (m − L)/m")
for name, d in dm.items():
    ax[0].axhline(d["chord_mean"], color="#c0392b" if "1m" in name else "#e67e22", lw=1,
                  label=f"damage x(T)/T, {name}: {d['chord_mean']:.4f} ± {d['chord_sem'] or 0:.4f}")
ax[0].axhline(0.252, color="grey", ls=":", lw=1, label="Wolfram 2019: ≈0.252 over first 10⁵ steps")
ax[0].set_xscale("log"); ax[0].set_ylim(0.235, 0.26); ax[0].set_xlabel("row m")
ax[0].set_ylabel("leftward speed (cells/step)"); ax[0].legend(fontsize=8); ax[0].set_title("Frozen-wedge frontier vs damage speed")
ax[1].plot(m[a::step], resid[::step], lw=0.4, color="#1f5fa8")
ax[1].set_xlabel("row m"); ax[1].set_ylabel("L(m) − linear fit (cells)"); ax[1].set_title("Frontier fluctuations about its line")
fig.tight_layout(); fig.savefig(f"{R}/frontier_vs_damage.png", dpi=130)
