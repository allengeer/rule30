#!/usr/bin/env python3
"""Charts for the 2026-10-27 X-Files (wedge, edge law, third dimension, long-range correlation)."""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

D = sys.argv[1] if len(sys.argv) > 1 else "/Volumes/BLIPS/rule30/2026-10-27_xfiles2"
INK, MUTED, GRID = "#1f2328", "#6e7781", "#d0d7de"
BLUE, ORANGE, GREEN, RED, PURPLE = "#2f6fdb", "#e07b39", "#2e9e6a", "#c9423a", "#7b5cc4"


def style(ax, title, xlabel, ylabel):
    ax.set_title(title, loc="left", fontsize=11, color=INK)
    ax.set_xlabel(xlabel, color=MUTED)
    ax.set_ylabel(ylabel, color=MUTED)
    ax.grid(True, color=GRID, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


fig, axes = plt.subplots(2, 2, figsize=(16, 11), dpi=130)
ax = axes.flat

# wedge boundary
ns, bound, frac = [], [], []
with open(os.path.join(D, "wedge/wedge.i64"), "rb") as f:
    while True:
        h = np.frombuffer(f.read(16), dtype=np.int64)
        if len(h) < 2:
            break
        n, Dn = h
        per = np.frombuffer(f.read(2 * Dn), dtype=np.int16)
        ns.append(n)
        bound.append(np.argmax(per == 0) / n)
        frac.append((per > 0).mean() * Dn / n)
a = next(ax)
a.plot(ns, bound, "o-", color=PURPLE, label="first aperiodic diagonal, d / n")
a.axhline(0.75, color=MUTED, ls="--", lw=0.8, label="3/4")
a.axhline(1.0, color=RED, ls=":", lw=1, label="centre column (d = n)")
a.set_ylim(0.7, 1.02)
a.legend(frameon=False)
style(a, "Frozen wedge · every diagonal left of ≈0.757n is periodic (periods 1–32)", "row n",
      "distance from left edge / n")

# edge law
r = np.fromfile(os.path.join(D, "edge_law_16m/rightrun.u16"), dtype=np.uint16).astype(int)
idx = np.arange(1, len(r) + 1)
v = np.array([(int(x) & -int(x)).bit_length() - 1 for x in 2 ** np.arange(0, int(np.log2(len(r))) + 1)])
rec = [r[2 ** k - 1] for k in range(int(np.log2(len(r))) + 1)]
a = next(ax)
a.plot(range(len(rec)), rec, "o-", color=ORANGE, label="rightmost run at row 2^v − 1 (A094604)")
for k, y in enumerate(rec):
    a.annotate(str(y), (k, y), textcoords="offset points", xytext=(-4, 6), fontsize=8, color=INK)
a.legend(frameon=False)
style(a, f"Edge law · run length depends only on v₂(n+1), verified for all {len(r):,} rows", "v = v₂(n+1)",
      "rightmost black run")

# third dimension
td = json.load(open(os.path.join(D, "thirddim/thirddim.json")))
a = next(ax)
labels, ints, weds, rands = [], [], [], []
for k in ("2", "3", "4"):
    labels.append(f"H({k}×{k} patch)")
    ints.append(td["T1_patch_entropy"][k]["interior"])
    weds.append(td["T1_patch_entropy"][k]["wedge"])
    rands.append(td["T1_patch_entropy"][k]["random_start"])
for k in ("2", "4", "8"):
    labels.append(f"H(coarse{k} | above)×10")
    ints.append(10 * td["T4_coarse_cond_entropy"][k]["interior"])
    weds.append(10 * td["T4_coarse_cond_entropy"][k]["wedge"])
    rands.append(10 * td["T4_coarse_cond_entropy"][k]["random_start"])
x = np.arange(len(labels))
a.bar(x - 0.27, ints, 0.27, color=BLUE, label="single-cell interior")
a.bar(x, rands, 0.27, color=MUTED, label="random-start Rule 30")
a.bar(x + 0.27, weds, 0.27, color=PURPLE, label="frozen wedge")
a.set_xticks(x, labels, fontsize=8)
a.legend(frameon=False, fontsize=9)
style(a, "Third dimension? Interior information = random-start to 4 decimals", "", "bits")

# variance-time
a = next(ax)
m = np.array([8, 1e3, 1e4, 1e5, 1e5, 1e6, 1e7, 1e8])
vr = np.array([0.9995, 0.9936, 1.0152, 0.9729, 1.0233, 1.047, 1.037, 0.844])
blocks = np.array([1.25e7, 1e5, 1e4, 1e3, 1e4, 1e3, 1e2, 10])
a.errorbar(m, vr, yerr=np.sqrt(2 / blocks), fmt="o", color=BLUE, capsize=3)
a.axhline(1, color=MUTED, ls="--", lw=0.8)
a.set_xscale("log")
style(a, "File 3 · centre column variance-time ratio (1 = no long-range correlation)", "block size (bits)",
      "Var(block sum) / (m/4)")

fig.suptitle("Rule 30 X-Files · 2026-10-27", x=0.01, ha="left", fontsize=15, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.97))
out = os.path.join(D, "xfiles_2026-10-27.png")
fig.savefig(out)
print(out)
