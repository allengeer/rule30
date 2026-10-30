# Rule 30 X-Files, 2026-10-23

Panel: Mulder (the anomaly is real), Scully (find the mundane cause), Data (quantify), Wesley (the unconventional angle), Picard (the ruling).

Instruments:
- `rule30.c`: 1M-row stats.
- `xfiles.c`: chain trace, random-start ring, edge scan.
- `crosscheck.py`: an independent numpy implementation.
- The Wolfram cloud kernel, plus the Wolfram Data Repository's 10⁹-bit centre column.

All runs are in this folder. Chart: `xfiles_2026-10-23.png`.

Coordinates:
- Row n starts at 0 with one black cell.
- Column s runs from −n (left edge) to +n (right edge).
- A chain is a maximal run of ≥ 2 adjacent 111 neighbourhoods, which is the seed of a pyramid (headstone).
- dl and dr are the distances from the left and right edges.

---

## X-File 1: The lumpy tail

**Hypothesis.** Scully: the excess of long chains is an edge effect. Mulder: it is structure in Rule 30 itself. Wesley: it is a pattern that drifts and grows.

**How the panel designed the test.** Scully wanted every long chain's position, not just a total count. Data asked for the right null: Rule 30 started from random rows rather than the "halve per step" formula. Picard required a second implementation before believing any count.

**Experiment.**
1. `xfiles trace 1000000` logged all 489,870 chains of length ≥ 18 with row, column and length. It also kept a histogram of every chain by length × log₂ distance to the nearer edge.
2. `xfiles random` ran four seeds × 20,000 steps on a 4,194,304-cell ring from fair random bits, 9.4 × 10⁹ chains in total.
3. Cross-checks:
   - Wolfram `CellularAutomaton` for the right-edge runs.
   - `crosscheck.py` (numpy) for rows 250,000–319,999.

**Findings.**
1. **Null confirmed.** In random-start Rule 30, each extra step of chain length halves the count, ratio 0.500 ± 0.01 out to length 26.
2. **Right edge: a power-of-two law.** All 3,906 long chains within 64 cells of the right edge sit on rows where n+1 is divisible by 2⁸ or more. The chain length is fixed by v, the exact power of 2 dividing n+1:
   - v = 8 → 22, 9 → 23, 10 → 25, 11 → 27, 12 → 32, 13 → 34
   - v = 14 → 35, 15 → 37, 16 → 39, 17 → 41, 18 → 46, 19 → 47
   - The counts are exactly what that rule predicts: 1,953, 977, 488, …, 1.
   - Wolfram confirms the rightmost run is 34 cells at row 4,095 (a chain of 32), and 1 cell at rows 4,094 and 4,096.
3. **Left wedge: a frozen region.** From the left edge to the line s ≈ −n/4, no chain of length 20 ever occurs, and none of length ≥ 22. Lengths 18, 19 and 21 are over-represented instead.
   - The boundary holds at s/n = −0.249, −0.249, −0.247 and −0.244 in rows 10k–100k, 100k–300k, 300k–600k and 600k–1M.
   - Those chains are persistent structures. Each sits at a fixed distance from the left edge and recurs **every 32 rows exactly**. Example: a chain of length 21 at dl = 125,351 recurs from row 167,280 onward.
   - 17 wedge diagonals carry at least one chain of length ≥ 18. Of these, **10 are persistent**: each holds 5,273–26,023 chains (6 of length 18, 2 of 19, 2 of 21), with every gap between recurrences a multiple of 32 rows. The other 7 carry exactly one chain each. Each persistent diagonal switches on when the −n/4 boundary sweeps past it (dl ≈ 0.75 n). (Recounted 2026-10-27.)
4. **Interior.** Between s = −n/4 and the right edge, every length ratio from 18 to 31 is ½ within noise, matching the random-start null.

**Conclusion.** Rule 30 contains no mysterious long-chain excess. The tail is three populations added together:
- The chaotic interior, which is exactly like a coin.
- A regular left wedge made of period-32 persistent structures.
- A right edge whose run length is set by the binary expansion of n+1.

Mulder gets a real structure, just not a supernatural one. Scully gets a mundane cause, and it is stranger than an "edge effect".

**Questions.**
- What sequence maps v to the edge chain length (22, 23, 25, 27, 32, 34, …)? Is it in OEIS? Can it be derived from the right edge's known nested structure?
- Is the wedge boundary speed exactly 1/4? Data notes the leftmost detections drift slightly (−0.249 → −0.244) and wants a direct measurement.
- Wesley: if the boundary ever reached s = 0, the centre column would become periodic. Problem 1 can be restated as "the frozen wedge never reaches the centre". Is that a known reformulation?
- Is the left wedge's period always 32, or 2ᵏ for different structures?
- Literature check: the regular left side of Rule 30 is well known (NKS). The −n/4 boundary and the period-32 structures may be too. Before the series claims novelty, search for "Rule 30 left boundary speed" and the Data Repository's persistent-structure sets.

---

## X-File 2: The missing 33

**Hypothesis.** Mulder: 33 is suppressed. Scully: it is one point inside File 1's pattern.

**How the panel designed the test.** Split the counts at 32–34 into populations and test the interior alone against the null.

**Experiment.** Population split from the File 1 trace.

**Findings.**
- The right-edge power-of-two law produces lengths 32 and 34 (122 and 61 chains) but never 33. Its sequence jumps 27 → 32 → 34.
- Interior only (dl, dr ≥ 64): lengths 30–35 = 40, 18, 10, **1**, 3, 0. The null predicts about 5 at length 33.
- P(X ≤ 1 | λ = 5) ≈ 0.04, and that is one test among about 15 tail lengths.

**Conclusion.** The "missing 33" is the right-edge sequence stepping over 33 while it inflates its neighbours. In the interior, 33 is unremarkable. Case closed. For the series, though: the edge sequence skips 33.

**Questions.** Why does the edge sequence jump 27 → 32 (a step of 5) and later 41 → 46? Are the steps related to the binary expansion of v?

---

## X-File 3: The centre column stays above ½

**Hypothesis.** Scully: a fair random walk stays on one side for long stretches, so this is expected. Mulder: the centre column has a persistent black bias.

**How the panel designed the test.** Data said our 1M rows could not settle it, so we predicted out of sample. If the walk is fair, the 10⁹-bit dataset should show it returning to a tie within about 1–2 × 10⁶ more steps with high probability. Picard asked for bit-level confirmation, since block checkpoints can miss a dip.

**Experiment.**
- In the Wolfram cloud kernel: the 10⁹-bit Data Repository dataset, counted in 10,000 blocks of 10⁵ bits.
- Bit-level scans of every block whose running excess started below 2,500, plus the first 1.1 × 10⁷ bits in full.
- The million-bit dataset and our engine agree exactly (500,768 ones at 10⁶).

**Findings.**
- The last tie (ones = n/2) is at **n = 195,112**, matching our engine.
- From there to n = 10⁹ the excess never returns to 0:
  - The smallest excess after 10⁶ is +349, at n = 3,428,516.
  - The smallest in the 286M region is +1,932.5. A block that started above 2,500 would need a drop of about 15σ.
- Final count: 500,025,038 ones, z = +1.58.
- Block counts: χ² = 10,232.7 on 10,000 df (z ≈ +1.6); lag-1 and lag-2 autocorrelation +0.008 and +0.009 (each < 1σ).
- Fair-walk probabilities:
  - Last tie before n = 195,112 out of 10⁹, by the arcsine law: (2/π)·arcsin√(1.95 × 10⁻⁴) ≈ **0.9%**.
  - The out-of-sample prediction: starting at +768 at 10⁶, never hitting 0 by 10⁹ has probability ≈ 2Φ(0.049) − 1 ≈ **3.9%**.

**Conclusion.** Scully's hypothesis predicted a return to a tie, and it didn't come in 999 million steps (p ≈ 0.04 on the prediction we wrote down first). Notable, but not decisive. The overall bias is z = 1.58, not significant, and the block statistics look fair. The file stays open.

**Questions.**
- Does the walk tie again beyond 10⁹? We'd need more bits; Wolfram's own records go further.
- Is the slight χ² excess (z ≈ 1.6) real long-range correlation? Test it with spectral and variance-time analysis over the 10⁹ bits.
- Does the wedge, or the period-32 structures, couple to the centre column at all? It sits ¼ n to the right of the boundary.

---

## X-File 4: What should the active-pyramid density be?

**Hypothesis.** Wesley: measure it instead of deriving it. If Rule 30's interior behaves like Rule 30 run from random bits, the single-cell density should equal the random-start density.

**How the panel designed the test.** Run the identical pyramid bookkeeping on random-start rings, after a warm-up, across four seeds.

**Experiment.** `xfiles random 4194304 20000 2000 seed` for seeds 1–4.

**Findings.**
- Random start: active/width = 0.078125 ± 0.000001 = **5/64 exactly**; new/width = 0.031250 = 1/32; density = 0.500000.
- The steady state is reached within about 5 steps: 0.0313, 0.0626, 0.0742, then 0.0780.
- Single cell, rows ≥ 500k: active = 0.078143, new = 0.031255. The small excess is plausibly the wedge's extra long chains.

**Conclusion.** The density matches the fair-coin model to within 0.02%. The 5/64 is clean enough that an exact derivation should exist. Closed.

**Questions.** Derive the 5/64. Data suggests it is the expected number of distinct pyramid runs per cell under the Bernoulli measure, which Rule 30 preserves. Does the left wedge alone have a different density?

---

## X-File 5: The edge bias

**Hypothesis.** Scully: the edges have fixed structure that decays within a few cells. Mulder: if the bias reaches the pyramids, it explains File 1.

**How the panel designed the test.** Measure E along each diagonal at fixed distance from each edge, out to 2¹⁸ cells, and find where the rms z falls to the fair-coin value of 1.

**Experiment.** `xfiles edge 1000000 262144 16`: 62,500 sampled rows, both edges.

**Findings.**
- **Right edge:** biased for the first ~60 cells, which hold the power-of-two structure. From 64 to 262,144 cells, rms z = 1.00 ± 0.04, a perfect coin.
- **Left edge:** rms z ≈ 250 at every distance out to about 65,000 cells, meaning each diagonal is essentially constant or periodic. It is still about 140 at 131k–262k cells, where the diagonals have only recently been swept into the wedge.

**Conclusion.** The two edges are opposites. The right edge is a thin crust of binary structure over a perfect coin. The left edge is the frozen wedge from File 1, growing at ¼ cell per row and covering 37.5% of the triangle. Scully's "decays within a few cells" is true on the right and false on the left. Mulder's link is confirmed: the left bias and the File 1 lumps are the same object.

**Questions.**
- Measure the wedge directly by period-testing each diagonal. Is everything left of −n/4 exactly periodic?
- What generates the period-32 structures at the moving boundary? It may be the left analogue of the right edge's 2-adic law.

---

## Instrument checks
- `xfiles.c` histograms equal `rule30.c` histograms exactly.
- Our 10⁶ centre bits equal the Wolfram Data Repository million-bit set.
- The last tie at n = 195,112 agrees between our engine and Wolfram.
- The right-edge runs agree with Wolfram `CellularAutomaton`.
- `crosscheck.py` against `chains.i64` for rows 250,000–319,999: both list the same 20,201 chains, with none missing on either side. That includes the 4,486 / 1,470 / 2,930 chains of length 19 / 20 / 21.
