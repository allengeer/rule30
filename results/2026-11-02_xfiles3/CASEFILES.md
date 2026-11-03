# Rule 30 X-Files, 2026-11-02

Panel: Mulder, Scully, Data, Wesley, Riker (chair). Captain Picard's standing orders remain in force: no count is believed without a second implementation, and no finding without replication.

Chart: `frontier_vs_damage.png`. Numbers: `summary.json`. Coordinates as before: row m, d = distance from the left edge, speeds are cells per step leftward in absolute coordinates.

---

## X-File 1c: Is the frozen wedge the chaos's shadow?

**Hypothesis.** On 2026-10-27 the wedge frontier (0.2433) and the left damage speed (0.2450 ± 0.0009) agreed only to about 2σ.
- Scully: they are the same speed, and the wedge is simply where the chaotic core's influence hasn't reached.
- Mulder: the gap is real, and something lives between the two lines.

**How the panel designed the test.**
- Data: remove both possible biases. Measure the frontier with no look-back window, and the damage speed with far more seeds and far longer runs.
- Riker (the Captain's rule): both new instruments are checked against an independent numpy implementation before any number is read.

**Experiment.**
- `xfiles frontier 10000000`: every wedge period divides 32, so diagonal d breaks at row m iff cell(m, d) ≠ cell(m−32, d). L(m) is the leftmost break: one exact frontier value per row, with no window. 49 minutes for 10⁷ rows.
- `xfiles damage T seed`: a bit-parallel C rewrite of `damage.py`. One cell is flipped at 0.3 W of a random ring of 2T+128 cells, and the leftmost difference is tracked. 64 seeds × 10⁶ steps and 256 seeds × 10⁵ steps.
- `frontier_check.py` recomputes both modes in numpy: 0 mismatches over 3,968 frontier rows and 2 × 5,000 damage steps. The check caught a bug in the first damage setup, where the speed-1 right-hand front wrapped around the ring and posed as the left front. The flip site was moved so it can't wrap. The 10-27 `damage.py` used walls rather than a ring and was not affected.

**Findings.**
1. **The window was never the problem.** At row 10⁶ the window-free frontier gives L/m = 0.75732, against 0.7572 from the windowed method.
2. **Damage speed is flat in time.** The mean front position x(t)/t is 0.2430–0.2434 at every t from 10⁴ to 10⁶:

   | run | seeds | x(T)/T | line fit, last ¾ |
   |---|---|---|---|
   | 10⁵ steps | 256 | 0.24326 ± 0.00019 | 0.24341 ± 0.00023 |
   | 10⁶ steps | 64 | **0.24294 ± 0.00013** | 0.24278 ± 0.00016 |

   So the 10-27 value of 0.2450 ± 0.0009 (64 × 4×10⁴) was about 2σ high: a fluctuation, not a trend.
3. **The frontier settles onto the damage line.**

   | m | 10⁵ | 2×10⁵ | 5×10⁵ | 10⁶ | 2×10⁶ | 5×10⁶ | 10⁷ |
   |---|---|---|---|---|---|---|---|
   | (m − L)/m | 0.2529 | 0.2481 | 0.2450 | 0.2427 | 0.2430 | 0.2429 | 0.2433 |

   Over rows 2.5×10⁶–10⁷: chord speed 0.24333, line fit 0.24356, 8 independent blocks 0.2434 ± 0.0004. Against the damage speed 0.24294 ± 0.00013, the difference is 0.0004–0.0006, about 1σ of the frontier's own scatter.
4. **The early excess belongs to the frontier alone.** Between 10⁴ and 10⁵ rows the frontier runs at about 0.255, which matches Wolfram's 2019 figure of "about 0.252 over the first 100,000 steps". The damage speed shows no such excess at any t. From about 7×10⁵ rows on, the frontier sits on the damage line.
5. **Fluctuations.** The damage front's spread grows as t^0.50 (sd 327 cells at 10⁵, 1,043 at 10⁶), so it wanders like a random walk. The frontier wanders ±1,700 cells about its line over 2.5×10⁶–10⁷ rows.

**Conclusion.** Scully's picture holds. Measured cleanly, the frozen-wedge frontier and the left damage speed agree at about 0.243, to within the frontier's own scatter. The wedge is the region the chaotic core's influence has not yet reached, a light-cone shadow. The 2σ gap was a low-statistics damage measurement.
- Correction, on the record: in the interim report the panel said the damage speed slows with run length and that both speeds shared an early slowdown. Both claims were wrong. The longer data show the damage speed flat in time; only the frontier had the early excess.
- Wolfram's 0.252 is real, but it is the frontier's early transient, not its asymptotic speed.

**Replication (Riker).** 4 fresh seeds × 10⁷ damage steps (`damage_10m/`, 100 minutes), the frontier's own depth:

| t | 10⁶ | 2×10⁶ | 5×10⁶ | 10⁷ |
|---|---|---|---|---|
| x(t)/t, mean of 4 | 0.24344 | 0.24350 | 0.24345 | **0.24342 ± 0.00022** |

- The damage speed is flat from 10⁶ to 10⁷, as it was from 10⁴ to 10⁶.
- At 10⁷ it equals the frontier's 0.2433–0.2434 to within 0.0001.
- These 4 seeds sit 0.0005 above the 64-seed 10⁶ mean at every t, including t = 10⁶. That is seed-to-seed scatter (about 1σ for a 4-seed mean), not a trend.
- Pooled, the left damage speed is 0.2430–0.2434, and the frontier sits inside that range. The finding replicates.

**Questions.**
- Wesley: what is the early transient? Between 10⁴ and 10⁶ rows the frontier sat up to 1,200 cells left of where the chaos could reach from the origin. Those diagonals are outside the light cone but not yet periodic, so they should be still-settling transients. Test against the transient lengths in OEIS A363346.
- Is there an exact value for the left damage speed (Rule 30's left Lyapunov exponent)? Read Gravner et al.'s Lyapunov-profile table and the Grassberger reference.
- Why does the frontier wander less than a damage front? Over 7.5×10⁶ rows it stays within ±1,700 cells of its line, while a single damage front spreads by about 1,000 cells in 10⁶ steps. Measure its fluctuation exponent directly.
- The centre sits about 0.243 n inside the chaos and the gap grows linearly, so this route does not threaten Problem 1.

---

## Instruments
- `xfiles frontier N dir` writes `frontier.i64`: L(m) for every row.
- `xfiles damage T seed dir` writes `damage_<seed>.i32`: front distance after each step.
- `frontier_check.py`: independent numpy check of both.
- `frontier_analysis.py`: speeds, block errors, `summary.json`, chart.
