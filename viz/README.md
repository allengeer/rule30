# Rule 30 Observatory

A local browser viewer for the Rule 30 triangle and the 10⁹-bit centre column.

```
./serve.sh            # http://127.0.0.1:8730/
```

Definitions match `../rule30.c` and `../xfiles.c`:
- Rows start at n = 0 with one black cell, and columns s run from −n to n.
- A chain is a maximal run of ≥ 2 adjacent 111 neighbourhoods.
- A chain's seed cells are its own cells plus the 110 cell to its right.
- Pyramid ("tombstone") membership persists through 000 neighbourhoods.

## Layers and the scale each reaches

| Layer | Data | Reach |
|---|---|---|
| Live JS engine (`static/engine.js`) | bit-packed 32-bit words, all rows in memory | 16,384 rows (10⁴ rows evolve in ~0.1 s) |
| Triangle, full resolution | on-demand C (`r30lib.c`), evolved from a checkpoint every 1,024 rows | anywhere in the 10⁶ triangle, ~0.3 s per 256² tile |
| Triangle, zoomed out | `level_<B>.bin` block grids, B = 128 … 32768 | the whole 10⁶ triangle |
| Long chains as dots | `chains.i32`: all 7,447,888 chains of length ≥ 14 | 10⁶ rows, at any zoom |
| Pattern search, triangle | on-demand C | up to 512 cells/px; coarser is refused |
| Centre column | Wolfram Data Repository bits, plus `lvl_k` min/max/end excess and per-2¹⁶ run histograms | all 10⁹ bits, zoomable to single bits |

A 10⁹-row triangle is 10¹⁸ cells and cannot be computed. Above 10⁶ rows only the centre column exists.

## Files

- `r30gen.c`: the tile generator. `r30gen N K B LMIN HB out` writes checkpoints, block levels, the chain index, the histograms and the per-row records.
- `r30lib.c`: the window renderer and the centre-column scans, loaded through ctypes.
- `centre_prep.c`: extracts and indexes the 10⁹ bits from the Data Repository MX. The raw bytes start at offset 239.
- `server.py`: static files plus a JSON/binary API (stdlib and numpy).
- `static/`: the UI (`app.js`), the engine (`engine.js`) and the styles.
- `validate.py` and `validate_js.mjs`: the checks listed below.

Data lives on BLIPS:
- `/Volumes/BLIPS/rule30/viz/tiles/1m`: 3.5 GB, built in 3m42s.
- `/Volumes/BLIPS/rule30/viz/centre`: 140 MB.
- `/Volumes/BLIPS/rule30/datarepo/billion.mx`: the 125 MB download.

## Rebuild

```
make
build/centre_prep /Volumes/BLIPS/rule30/datarepo/billion.mx 239 1000000000 /Volumes/BLIPS/rule30/viz/centre
build/r30gen 1000000 1024 128 14 1024 /Volumes/BLIPS/rule30/viz/tiles/1m
../.venv/bin/python validate.py 1m
node validate_js.mjs 10000 http://127.0.0.1:8730 test20k   # with the server running
```

## Checks (all pass, 2026-10-27)

- `r30gen` agrees with `rule30.c`'s `rows.i64` on centre bit, ones, new chains and active pyramids for all 10⁶ rows. Its chain-length histogram equals `hist.i64` (31,255,819,693 chains).
- `r30gen` chains of length ≥ 18 equal `xfiles.c`'s `chains.i64` exactly: 489,870 chains.
- `r30_window` agrees with a plain numpy evolution on cells, chain cells, pyramid labels, `11011` matches and histograms. That covers 12 windows at 1 cell/px and the block counts at 2, 8 and 64 cells/px, including windows that straddle the edges.
- On-demand tiles, summed 2×2, equal the precomputed `level_128`. Level totals equal the per-row records, and every coarser level preserves totals.
- The JS engine agrees with `rows.i64` (10⁴ rows) and `chains.i64`. Its tiles match the C server's byte-for-byte on 42 tiles, with and without patterns.
- The Data Repository bits equal `rule30.c`'s centre column for the first 10⁶ bits. The totals are 500,025,038 ones and 339 ties, the last at 195,112. All 1,000 entries of `walk_every_1e6.csv` match.

Pyramid labels in a window are exact for chains shorter than about 250. The renderer starts 128 rows above the window, and a pyramid shrinks by at least one cell per side per row.
