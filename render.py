#!/usr/bin/env python3
"""Render a dense sim.py run: the Rule 30 pattern, E(cell) and pyramid-count heatmaps."""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

Image.MAX_IMAGE_PIXELS = None
RULE30 = np.array([(30 >> k) & 1 for k in range(8)], dtype=np.uint8)


def pattern(N):
    W = 2 * N + 1
    img = np.zeros((N, W), dtype=np.uint8)
    row = np.zeros(W + 2, dtype=np.uint8)
    row[N + 1] = 1
    for n in range(N):
        img[n] = row[1:-1]
        nb = (row[:-2] << 2) | (row[1:-1] << 1) | row[2:]
        row[1:-1] = RULE30[nb]
    return img


def sampled(path, dtype, N, rows_out, cols_out, fill, zscore=False):
    """Downsample the triangular per-row file to rows_out x cols_out by row sampling and column-bin means."""
    data = np.memmap(path, dtype=dtype, mode="r")
    W = 2 * N - 1
    edges = np.linspace(0, W, cols_out + 1).astype(int)
    out = np.full((rows_out, cols_out), fill, dtype=np.float64)
    for j, n in enumerate(np.linspace(0, N - 1, rows_out).astype(int)):
        full = np.full(W, np.nan)
        vals = data[n * n:n * n + 2 * n + 1].astype(np.float64)
        if zscore:  # (E - 1/2) / sqrt(1/(4m)), m = rows the column has been in the cone
            m = n - np.abs(np.arange(-n, n + 1)) + 1
            vals = (vals - 0.5) * 2 * np.sqrt(m)
        full[N - 1 - n:N + n] = vals
        with np.errstate(invalid="ignore"):
            sums = np.add.reduceat(np.nan_to_num(full), edges[:-1])
            cnt = np.add.reduceat(~np.isnan(full), edges[:-1])
            out[j] = np.where(cnt > 0, sums / np.maximum(cnt, 1), np.nan)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="out_1gb")
    ap.add_argument("--img", default="img")
    args = ap.parse_args()
    os.makedirs(args.img, exist_ok=True)
    N = int(os.path.getsize(os.path.join(args.run, "pyr.u16")) // 2) ** 0.5
    N = int(round(N))

    pat = pattern(N)
    Image.fromarray((1 - pat) * 255).convert("1").save(os.path.join(args.img, "rule30_full.png"), optimize=True)
    # preview: 4x4 block mean as grey
    k = max(1, pat.shape[1] // 3200)
    h, w = pat.shape[0] // k * k, pat.shape[1] // k * k
    prev = pat[:h, :w].reshape(h // k, k, w // k, k).mean(axis=(1, 3))
    Image.fromarray((255 * (1 - prev)).astype(np.uint8)).save(os.path.join(args.img, "rule30_preview.png"))
    centre = pat[:, N]
    del pat

    E = sampled(os.path.join(args.run, "e.f32"), np.float32, N, 1200, 2400, np.nan, zscore=True)
    P = sampled(os.path.join(args.run, "pyr.u16"), np.uint16, N, 1200, 2400, np.nan)
    ext = [-(N - 1), N - 1, N - 1, 0]
    for data, name, kw, label in [
        (E, "E_heatmap.png", dict(cmap="RdBu_r", vmin=-1.5, vmax=1.5), "mean z over ~11 columns, z = (E - 1/2) * 2 sqrt(m)"),
        (np.log10(1 + P), "pyramids_heatmap.png", dict(cmap="magma"), "log10(1 + pyramids covering column)"),
    ]:
        fig, ax = plt.subplots(figsize=(16, 8), dpi=150)
        im = ax.imshow(data, aspect="auto", extent=ext, interpolation="nearest", **kw)
        ax.set_xlabel("column offset from centre")
        ax.set_ylabel("row")
        fig.colorbar(im, ax=ax, label=label)
        fig.tight_layout()
        fig.savefig(os.path.join(args.img, name))
        plt.close(fig)
    print(f"rendered {N} rows; centre ones = {int(centre.sum())}")


if __name__ == "__main__":
    main()
