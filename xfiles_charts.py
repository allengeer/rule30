#!/usr/bin/env python3
"""Charts for the 2026-10-23 X-File investigations (reads the BLIPS run directory)."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

D = sys.argv[1] if len(sys.argv) > 1 else "/Volumes/BLIPS/rule30/2026-10-23_xfiles"
INK, MUTED, GRID = "#1f2328", "#6e7781", "#d0d7de"
BLUE, ORANGE, GREEN, RED, PURPLE = "#2f6fdb", "#e07b39", "#2e9e6a", "#c9423a", "#7b5cc4"


def style(ax, title, xlabel, ylabel):
    ax.set_title(title, loc="left", fontsize=11, color=INK)
    ax.set_xlabel(xlabel, color=MUTED)
    ax.set_ylabel(ylabel, color=MUTED)
    ax.grid(True, color=GRID, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


ch = np.fromfile(os.path.join(D, "trace_1m/chains.i64"), dtype=np.int64).reshape(-1, 3)
n, s, L = ch.T
v = s / np.maximum(n, 1)
right = (n - (s + L - 1)) <= 64
wedge = (v < -0.25) & ~right
rest = ~wedge & ~right

fig, axes = plt.subplots(2, 2, figsize=(16, 11), dpi=130)
ax = axes.flat

# X-File 1/2: where the tail comes from
a = next(ax)
k = np.arange(18, 50)
rh = sum(np.fromfile(os.path.join(D, f"random/seed{i}/random_hist.i64"), dtype=np.int64) for i in range(1, 5))
for mask, colour, label in [(rest, BLUE, "interior, right of the −n/4 line"),
                            (wedge, PURPLE, "left wedge (s/n < −1/4)"),
                            (right, ORANGE, "right edge (≤ 64 cells)")]:
    c = np.bincount(L[mask], minlength=60)[k]
    sel = c > 0
    a.semilogy(k[sel], c[sel], "o-", color=colour, ms=4, lw=1, label=label)
base = np.bincount(L[rest], minlength=60)[18]
a.semilogy(k, base * 0.5 ** (k - 18), color=MUTED, ls="--", label="random-start Rule 30 (measured ×½ per step)")
a.axvline(33, color=RED, lw=0.8, ls=":")
a.text(33.3, 2e5, "33", color=RED)
a.legend(frameon=False, fontsize=9)
style(a, "X-Files 1–2 · chain lengths ≥ 18 by where they occur (1M rows)", "chain length (adjacent 111s)", "count")

# the wedge
a = next(ax)
for length, colour in [(19, GREEN), (20, RED), (21, PURPLE)]:
    m = (L == length) & ~right
    idx = np.flatnonzero(m)[:: max(1, m.sum() // 6000)]
    a.scatter(v[idx], n[idx], s=1, color=colour, alpha=0.5, label=f"length {length}")
a.axvline(-0.25, color=INK, lw=1, ls="--")
a.text(-0.245, 50000, " s = −n/4", color=INK)
a.set_xlim(-1, 1)
a.legend(frameon=False, markerscale=8, loc="upper right")
style(a, "X-File 1 · no length-20 chain ever appears left of s = −n/4", "position s / n (−1 = left edge)", "row n")

# X-File 3: the billion-bit walk
a = next(ax)
S = np.loadtxt(os.path.join(D, "xfile3_billion/walk_every_1e6.csv"), delimiter=",")
m = np.arange(1, len(S) + 1) * 1e6
z = S / (np.sqrt(m) / 2)
a.plot(m, z, color=BLUE, lw=1)
for kk in (-2, 0, 2):
    a.axhline(kk, color=RED if kk else MUTED, lw=0.8, ls="--")
a.set_xscale("log")
a.set_ylim(-3, 3.5)
style(a, "X-File 3 · centre column, 10⁹ bits: last tie at n = 195,112, then never again",
      "n (centre-column cells)", "z = (ones − n/2)/(√n/2)")

# X-File 5: edge reach
a = next(ax)
e = np.fromfile(os.path.join(D, "edge_wide/edge_wide.i64"), dtype=np.int64).reshape(3, -1)
d = np.arange(e.shape[1])
edges = np.unique(np.geomspace(1, len(d), 60).astype(int))
for side, colour, label in [(0, PURPLE, "left edge"), (1, ORANGE, "right edge")]:
    zz = ((e[side] / e[2]) - 0.5) * 2 * np.sqrt(e[2])
    rms = [np.sqrt(np.mean(zz[a0:a1] ** 2)) for a0, a1 in zip(edges[:-1], edges[1:])]
    a.loglog(edges[:-1], rms, "o-", ms=3, color=colour, label=label)
a.axhline(1, color=MUTED, ls="--", lw=0.8, label="fair coin (rms z = 1)")
a.legend(frameon=False)
style(a, "X-File 5 · how far each edge's bias reaches (rows sampled every 16th)",
      "distance from edge (cells)", "rms z of diagonal E")

fig.suptitle("Rule 30 X-Files · 2026-10-23", x=0.01, ha="left", fontsize=15, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.97))
out = os.path.join(D, "xfiles_2026-10-23.png")
fig.savefig(out)
print(out)
