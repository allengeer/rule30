# Rule 30 X-Files, 2026-10-27

Panel: Mulder, Scully, Data, Wesley, Picard. Data is building the visual explorer in a separate session (`rule30-viz`).

Chart: `xfiles_2026-10-27.png`. Coordinates as in the 2026-10-23 files: row n, column s ∈ [−n, n], and d = distance from the left edge.

---

## X-File 1b: The frozen wedge, mapped

**Hypothesis.**
- Wesley: everything left of about −n/4 is exactly periodic, and the boundary moves at exactly ¼.
- Scully: it's the known regular left side of Rule 30 (NKS), so check the literature first.
- Mulder: if the boundary moves, which way?

**How the panel designed the test.** Data: don't infer periodicity from chains, test every diagonal directly. At each checkpoint, take every diagonal's last 256 cells and find its smallest period (up to 128). The boundary is the first diagonal with no period. Picard: and search the literature in parallel.

**Experiment.** `xfiles wedge 1000001 100000 256` checked 10 checkpoints, about 1.1 × 10⁷ diagonals in total. Literature: web search, OEIS A094605 / A363345 / A363346, and Rowland (2006).

**Findings.**
- At every checkpoint, **every** diagonal from the left edge to the boundary is periodic, and every diagonal beyond about 0.80 n is not.
- The periods are powers of two: 1, 2, 4, 8, 16, 32 (diagonal counts at row 10⁶: 16, 10, 56, 668, 87,133, 669,279). The period-16 band stops growing near 87k diagonals, while the period-32 band grows linearly. No period 64 has appeared yet.
- The boundary (first aperiodic diagonal ÷ n) **creeps upward**:

  | n | 100k | 200k | 300k | 400k | 500k | 600k | 700k | 800k | 900k | 1M |
  |---|---|---|---|---|---|---|---|---|---|---|
  | boundary | .7454 | .7511 | .7510 | .7536 | .7547 | .7554 | .7567 | .7567 | .7570 | .7572 |

  That is about +0.006 per tenfold increase in n after 200k. The 256-row window could bias this by at most about 0.002 at 100k.
- **Literature:**
  - The left-diagonal periodicity and period doubling are known: Wolfram's NKS, proved by Rowland, *Local Nested Structure in Rule 30*, Complex Systems 16 (2006).
  - The eventual periods and transients of the left diagonals are OEIS A363345 and A363346.
  - We did not find a published boundary speed.

**Conclusion.**
- The wedge is the known periodic left side. Our new contribution is a measurement: its frontier sits at d ≈ 0.757 n at row 10⁶ and is still creeping toward the centre.
- Mulder's question has an answer: the wedge grows *faster* than n.
- At the observed rate it would need roughly 10³⁰ rows to reach the centre (d = n). It's a slow creep, but it bears directly on Problem 1: the centre column stays aperiodic exactly as long as this frontier never reaches it.

**Questions.**
- Is the creep logarithmic and unbounded, or does it settle to a constant just above ¾?
- Extend to 10⁷ rows. Stepping 10⁶ rows takes about 3.3 s and cost grows as n², so 10⁷ rows is about 6 minutes plus the periodicity checks. Running: `wedge_10m/`.
- When does period 64 first appear, and how does that fit the A094605-style doubling?
- Can our boundary data be compared against A363346 (transient lengths)?

### Update: out to 10⁷ rows (`wedge_10m/`)

**Findings.**
- The creep **stopped**. The frontier d/n was 0.7572 at 1M, then held flat: 0.7569, 0.7572, 0.7570, 0.7570, 0.7569, 0.7567, 0.7567, 0.7567 and 0.7567 at 2M through 10M.
- At this point the frontier is a straight line, moving at constant speed: s/n ≈ −0.2433.
- Still no period-64 diagonal at 10⁷.
- The period-32 band grows linearly, about 757k diagonals per 10⁶ rows. The period-16 band keeps creeping up slowly (87,133 → 87,257).

**The panel's next test.**
- **Data:** a constant frontier speed near 0.24 resembles Rule 30's leftward damage-spreading speed.
- **Scully:** if they are the same, the wedge needs no mystery. Its frontier is simply as far left as the chaotic core can send any influence. So measure the damage speed independently.

**Experiment.** `damage.py`: random row, flip one cell, run both copies, and track the leftmost differing cell. A linear fit over the last ¾ of each run.
- 16 seeds × 20,000 steps: speed 0.2445 ± 0.0022.
- 64 seeds × 40,000 steps: speed **0.2450 ± 0.0009**.

**Conclusion.** The frozen frontier (0.2433) and the left damage speed (0.2450 ± 0.0009) agree within about 2σ. The likely explanation is that the wedge is the region the chaotic core's influence hasn't reached, which makes it a light-cone shadow rather than an independent structure. Mulder's creep was the frontier settling into its asymptotic line, and it is not heading for the centre.
- For Problem 1 this cuts the other way from what Mulder hoped. The centre sits deep inside the chaotic core, roughly a quarter of n from the frontier, and the gap widens linearly.

**Questions.**
- Is the 2σ gap real? Measure the frontier with a method free of window bias, and the damage speed with more seeds.
- Is the damage speed known exactly, as Rule 30's left Lyapunov exponent? Literature check.

---

## X-File 2b: The edge law

**Hypothesis.**
- Scully: the right-edge power-of-two law is already known.
- Data: then check whether it holds everywhere, not just on long chains.

**How the panel designed the test.** Compute the rightmost black run for every row, not just rows with chains of 18 or more, and test the exact claim: the run length depends only on v₂(n+1).

**Experiment.**
- `xfiles rightrun` for 2²⁰ rows (3.3 s), then 2²⁴ rows.
- OEIS A094603 / A094604.

**Findings.**
- It is **OEIS A094603** (rightmost black run per row), with records **A094604** at rows 2ᵐ. Rowland (2006) explains why: at row 2ⁿ the initial condition reappears on the right and the automaton "begins again" locally.
- A094603 carries an open conjecture (Gómez Calderón, 2023): a(n) = A094604(v₂(n)).
- **Verified for every n ≤ 2²⁴ = 16,777,216** (the OEIS b-file has 10⁴ terms). In `edge_law_16m/rightrun.u16` (33,554,432 bytes, complete), every row with v₂(n+1) = v has a single run length, for v = 0…24. The records for v = 21…24 are 54, 55, 58, 60, matching A094604.
- Our records match A094604's published terms exactly: 1, 3, 4, 6, 7, 9, 15, 16, 24, 25, 27, 29, 34, 36, 37, 39, 41, 43, 48, 49, 51.

**Conclusion.** Not a discovery. It's a known structure, and we've confirmed an open OEIS conjecture over far more terms. Those are computational checks, not a proof; a proof probably follows from Rowland's reversibility argument.

**Questions.**
- Does Rowland's theorem already imply the conjecture? Read the paper and draft the argument (in Lean, per the Prize skill, if it's worth it).
- Is a submission to the OEIS b-file worth doing?

---

## X-File 6: Is the 2-D picture encoding a third dimension?

**Hypothesis.**
- Mulder: the triangle is a projection, and a hidden coordinate is organising it.
- Wesley: that coordinate is *scale*, the powers of two, given 2ᵏ periods on the left, v₂(n+1) on the right, and restarts at 2ⁿ.
- Scully: Rule 30 is deterministic from one row, so the second dimension is just time. Any "extra" dimension has to show up as measurable information or correlation that a plain 1-D rule wouldn't produce.

**How the panel designed the test.** Data listed four fingerprints a hidden coordinate would leave. Each is compared against the same rule run from random bits, which has no special origin and no hidden structure, and split into interior and wedge:
1. **Patch information:** the entropy of k×k patches. A hidden constraint lowers it.
2. **Row correlations:** corr(x_s, x_{s+Δ}) for Δ = 1…2048. Hidden coupling means non-zero correlation.
3. **Scale index:** do interior statistics depend on v₂(n+1)? That would mean scale acts as a coordinate.
4. **Renormalisation:** coarse-grain into k×k majority blocks (k = 2, 4, 8) and measure how predictable each coarse cell is from the coarse row above. A hidden structure across scales means lower entropy.

Picard: anything found must be replicated on a fresh sample before we believe it.

**Experiment.** `thirddim.py`: single-cell Rule 30 to 16,384 rows, plus a random ring of width 32,768, measured on rows 8,192–16,383. Then a replication of test 3 at 32,768 rows against 8 random-start seeds.

**Findings.**
1. **Patches.** Interior and random-start agree to the 4th decimal:
   - 2×2: 3.50015 vs 3.50001 bits.
   - 3×3: 5.82796 vs 5.82782.
   - 4×4: 8.12301 vs 8.12281.

   The wedge is slightly lower (8.10385 at 4×4).
2. **Correlations.** Interior rms correlation 0.00077, max 0.0025. The random-start rms of 0.00050 scales to the 0.00079 expected for the interior's smaller sample, so it's pure noise. The wedge is strongly correlated (max 0.039), being periodic.
3. **Scale.** At 16k rows, interior density by v₂(n+1) gave naive χ² = 16.7 (p ≈ 0.02). **It did not replicate:** at 32k rows χ² = 1.4. And the same statistic in random-start runs ranges from 4.5 to 28 (rows are correlated), so the first hint was noise.
4. **Renormalisation.** Conditional entropy of coarse cells, interior vs random:

   | k | interior | random |
   |---|---|---|
   | 2 | 0.8365 | 0.8370 |
   | 4 | 0.9415 | 0.9414 |
   | 8 | 0.9885 | 0.9884 |

   The wedge is more predictable at k = 4 (0.9233).

**Conclusion.** In the chaotic interior there is no detectable hidden dimension. By all four fingerprints it is indistinguishable from Rule 30 run on random bits, which is flat 1+1-D chaos. Wesley's scale coordinate is real, but **only at the boundaries**: the right edge (indexed by v₂(n+1)) and the frozen wedge (2ᵏ periods). Rule 30 is a fractal frame around a chaotic core. Mulder gets a hidden binary scale axis; Scully gets that it lives only in the frame.

**Questions.**
- Run the fingerprints closer to the wedge frontier. Does structure leak across the boundary (a "membrane")?
- Is the interior/frontier transition sharp, and does its width grow?

---

## X-File 3 (continued): Is the persistence correlation or luck?

**Hypothesis.** Mulder: black stays ahead because of a long-range memory. Scully: it's luck (about a 4% event).

**How the panel designed the test.** Long-range memory inflates the variance of block sums at large block sizes, giving a Hurst exponent above ½. Measure Var(sum)/(m/4) for m from 8 to 10⁸.

**Experiment.**
- Wolfram cloud kernel on the first 10⁸ bits of the Data Repository set (m = 8, 10³, 10⁴, 10⁵).
- The full 10⁹ walk for 10⁵, 10⁶, 10⁷ and 10⁸.

**Findings.** Variance ratios, with the fair-coin standard error in brackets:

| m | 8 | 10³ | 10⁴ | 10⁵ (first 10⁸) | 10⁵ (all) | 10⁶ | 10⁷ | 10⁸ |
|---|---|---|---|---|---|---|---|---|
| ratio | 0.9995 | 0.994 (±0.004) | 1.015 (±0.014) | 0.973 (±0.045) | 1.023 (±0.014) | 1.047 (±0.045) | 1.04 (±0.14) | 0.84 (±0.45) |

All are within 1.6σ of 1. There's no trend with m.

**Conclusion.** No long-range correlation. The persistence is best explained as an approximately 4% fluctuation. The file stays open only for what comes after 10⁹.

---

## Instruments
- `xfiles wedge` and `xfiles rightrun` are new modes. Bit-parallel stepping runs at about 3.3 s for 2²⁰ rows.
- `thirddim.py` is pure numpy, independent of the C engine.
- The Wolfram cloud kernel handled the billion-bit tests.
- All outputs are under this folder.
