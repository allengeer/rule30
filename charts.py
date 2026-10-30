#!/usr/bin/env python3
"""Chart the per-row statistics of a deep rule30.c run against the row index."""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FIELDS = ["row", "cbit", "cones", "cperiod", "cpyr", "cinpyr", "newpyr", "active", "total",
          "ones", "onesl", "onesr"] + [f"t{k:03b}" for k in range(8)]
INK, MUTED, GRID = "#1f2328", "#6e7781", "#d0d7de"
BLUE, ORANGE, GREEN, RED = "#2f6fdb", "#e07b39", "#2e9e6a", "#c9423a"


def load(run):
    r = np.fromfile(os.path.join(run, "rows.i64"), dtype=np.int64).reshape(-1, len(FIELDS))
    d = {f: r[:, i] for i, f in enumerate(FIELDS)}
    n = d["row"]
    period = np.full(len(n), -1, dtype=np.int64)
    for path in glob.glob(os.path.join(run, "period_*_of_*.i64")):
        p = np.fromfile(path, dtype=np.int64).reshape(-1, 2)
        p = p[p[:, 0] < len(n)]
        period[p[:, 0]] = p[:, 1]
    d["rperiod"] = period
    return d


def logidx(n, k=4000):
    return np.unique(np.geomspace(1, n - 1, k).astype(int))


def style(ax, title, xlabel="row n", ylabel=None, logx=True):
    ax.set_title(title, loc="left", fontsize=11, color=INK)
    if logx:
        ax.set_xscale("log")
    ax.set_xlabel(xlabel, color=MUTED)
    if ylabel:
        ax.set_ylabel(ylabel, color=MUTED)
    ax.grid(True, color=GRID, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="out_1m")
    ap.add_argument("--img", default="img")
    args = ap.parse_args()
    os.makedirs(args.img, exist_ok=True)
    d = load(args.run)
    n = d["row"].astype(np.float64)
    N = len(n)
    width = 2 * n + 1
    i = logidx(N)

    fig, axes = plt.subplots(4, 2, figsize=(16, 20), dpi=130)
    ax = axes.flat

    # 1. centre E with the fair-coin 95% band
    m = n + 1
    E = d["cones"] / m
    a = next(ax)
    a.fill_between(m[i], 0.5 - 1.96 * 0.5 / np.sqrt(m[i]), 0.5 + 1.96 * 0.5 / np.sqrt(m[i]),
                   color=BLUE, alpha=0.12, lw=0, label="fair coin, 95%")
    a.plot(m[i], E[i], color=BLUE, lw=1, label="centre column E")
    a.axhline(0.5, color=MUTED, lw=0.8, ls="--")
    a.set_ylim(0.45, 0.55)
    a.legend(frameon=False)
    style(a, "Problem 2 · centre-column E(black) so far", ylabel="E")

    # 2. centre deviation as a z-score (random walk scaled by sqrt n)
    z = (d["cones"] - m / 2) / (np.sqrt(m) / 2)
    a = next(ax)
    a.plot(m, z, color=BLUE, lw=0.6)
    for k in (-2, 2):
        a.axhline(k, color=RED, lw=0.8, ls="--")
    style(a, f"Problem 2 · z = (ones − n/2)/(√n/2)   final z = {z[-1]:+.2f}", ylabel="z", logx=False)

    # 3. periodicity
    a = next(ax)
    a.plot(m[i], d["cperiod"][i] / m[i], color=BLUE, lw=1, label="centre prefix period / (n+1)")
    ok = d["rperiod"] > 0
    j = i[ok[i]]
    a.plot(m[j], d["rperiod"][j] / width[j], color=ORANGE, lw=1, label="row period / (2n+1)")
    a.set_ylim(0, 1.1)
    a.legend(frameon=False, loc="lower right")
    style(a, "Problem 1 · smallest period relative to length (1 = no repetition)", ylabel="ratio")

    # 4. row densities
    a = next(ax)
    h = n > 0
    a.plot(n[h][::50], (d["onesl"][h] / n[h])[::50], color=BLUE, lw=0.5, label="left half")
    a.plot(n[h][::50], (d["onesr"][h] / n[h])[::50], color=ORANGE, lw=0.5, label="right half")
    a.axhline(0.5, color=MUTED, lw=0.8, ls="--")
    a.set_ylim(0.45, 0.6)
    a.legend(frameon=False)
    style(a, "Row density of black cells (every 50th row)", ylabel="fraction black", logx=False)

    # 5. neighbourhood triples
    a = next(ax)
    cmap = plt.get_cmap("tab10")
    cw = np.cumsum(width)
    for k in range(8):
        f = np.cumsum(d[f"t{k:03b}"]) / cw
        a.plot(n[i], f[i], lw=1, color=cmap(k), label=f"{k:03b}")
    a.axhline(1 / 8, color=MUTED, lw=0.8, ls="--")
    a.set_ylim(0.10, 0.15)
    a.legend(frameon=False, ncol=4, fontsize=8)
    style(a, "Neighbourhood (triple) frequencies, all rows up to n  (uniform = 1/8)", ylabel="fraction of cells")

    # 6. pyramids
    a = next(ax)
    a.plot(n[i], d["newpyr"][i] / width[i], color=BLUE, lw=1, label="new pyramids / width")
    a.plot(n[i], d["active"][i] / width[i], color=ORANGE, lw=1, label="active pyramids / width")
    a.plot(n[i], d["cpyr"][i] / m[i], color=GREEN, lw=1, label="centre pyramids / (n+1)")
    a.axhline(1 / 32, color=BLUE, lw=0.8, ls="--", label="random row: 1/32 new")
    a.axhline(1 / 8, color=GREEN, lw=0.8, ls="--", label="random row: 1/8 covered")
    a.legend(frameon=False)
    style(a, "Pyramids (chains of ≥2 adjacent 111s)", ylabel="rate")

    # 7. chain-length histogram vs random-row expectation
    a = next(ax)
    hist = np.fromfile(os.path.join(args.run, "hist.i64"), dtype=np.int64)
    k = np.arange(len(hist))
    sel = (k >= 2) & (hist > 0)
    a.semilogy(k[sel], hist[sel], "o", color=BLUE, ms=4, label="Rule 30")
    a.semilogy(k[sel], hist[2] * 0.5 ** (k[sel] - 2), color=MUTED, ls="--", label="random row (×½ per step)")
    a.legend(frameon=False)
    style(a, f"Pyramid chain-length distribution  ({hist.sum():,} pyramids)", xlabel="chain length (111s)",
          ylabel="count", logx=False)

    # 8. final-row column E as z-scores
    a = next(ax)
    snaps = sorted(glob.glob(os.path.join(args.run, "snap_*.bin")), key=lambda p: int(p.split("_")[-1][:-4]))
    last = snaps[-1]
    sn = int(last.split("_")[-1][:-4])
    raw = np.fromfile(last, dtype=np.uint8)
    w = 2 * sn + 1
    Es = raw[:4 * w].view(np.float32).astype(np.float64)
    c = np.arange(-sn, sn + 1)
    mm = sn - np.abs(c) + 1
    zc = (Es - 0.5) * 2 * np.sqrt(mm)
    bins = 2000
    edges = np.linspace(0, w, bins + 1).astype(int)
    zb = np.add.reduceat(zc, edges[:-1]) / np.diff(edges)
    a.plot(c[edges[:-1]], zb, color=BLUE, lw=0.6)
    for kk in (-2, 2):
        a.axhline(kk, color=RED, lw=0.8, ls="--")
    a.set_ylim(-6, 6)
    style(a, f"Per-column E z-score at row {sn:,} (binned means of ~{w // bins} columns)",
          xlabel="column offset from centre", ylabel="mean z", logx=False)

    fig.suptitle(f"Rule 30 · {N:,} rows", x=0.01, ha="left", fontsize=15, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    out = os.path.join(args.img, "rule30_1m_stats.png")
    fig.savefig(out)
    print(out)
    print(f"centre E={E[-1]:.6f} z={z[-1]:+.3f}; centre period={d['cperiod'][-1]} of {N}; "
          f"min row-period ratio (n>100)={np.min((d['rperiod'] / width)[100:][ok[100:]]):.4f}; "
          f"pyramids={d['total'][-1]:,}")


if __name__ == "__main__":
    main()
