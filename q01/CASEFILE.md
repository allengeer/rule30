# X-File 8: Period two (Q01), 2026-11-03

Panel: Mulder, Scully, Data, Wesley, Riker (chair). Captain Picard's standing orders in force.

## Where Problem 1 stands
- Wolfram Prize Problem 1 asks for a **proof** that the centre column never becomes periodic. No amount of computation wins it.
- Proven so far: Jen (1986), at most one column of a finite orbit is eventually periodic. Condrey ([arXiv:2609.09431](https://arxiv.org/abs/2609.09431), Sep 2026), no column is eventually constant (period 1); he names period 2 as the next open case.
- The [Rule 30 Prize Kit](https://github.com/WolframInstitute/Rule30PrizeKit) (Lean 4) reduces period 2 to **Q01** (`Rule30ZeroBoundary.Q01`) and proves it for starts of width ≤ 52. Integer form (`Q01_iff_integer`): with F(Y) = (Y>>1) ⊕ (Y ∨ (Y<<1)) and G = F∘F, every natural Y reaches an even number under some iterate of G.
- Q01 ⇒ the centre column is not eventually of period 2 (`no_eventual_period_two`).

## Experiment
- `escape.c`: brute force over every odd Y < 2⁴⁰, counting G-steps to an even number (22 min on 4 cores).
- `tree.c`: bit-tree search. Bit 0 after s steps depends only on bits 0..2s, so bits are fixed from the bottom and a branch is dropped once it reaches an even number. 256-bit arithmetic with masking, so nothing overflows. Bit lengths 41–66, step cap W + 20 (94 min on 4 cores).
- Cross-checks (Captain's rule):
  - `tree.c` reproduces `escape.c`'s record and holder count at W = 11, 20, 30, 38, 40.
  - Its leaf count at even W is exactly 2^(W/2 − 1), consistent with the Kit's proved `count_constrained` (2^⌊W/2⌋ constrained starts of width ≤ W).

## Findings
1. **No counterexample below 2⁶⁶.** Every odd Y of bit length ≤ 66 reaches an even number. That is Q01 for all starts of width ≤ 66, against the Kit's Lean-proved 52. This is a computation, not a kernel-checked proof.
2. **The escape time is linear in the bit length**, with record holders listed in `escape_records_2e40.txt` and `tree_41_66.txt`:
   - Never more than W + 4 steps; that maximum is reached only at W = 11 (Y = 0x5d5).
   - From W = 41 to 66 the records sit between W − 9 and W − 2, so the slope is below 1: the record is 60 steps at W = 65.
3. **Record holders end in the stripe.** Many records have low bits …5555 (0x35555, 0x19b555, 0xd81ba03555555, 0x1e93fab5a49235555). The stripe …10101 is the ZBL's fixed point (`zrow_stripe`) and never reaches an even number. The slowest finite starts are long stripe prefixes with a head that delays the stripe's collapse.

## Conjecture E (effective Q01)
Every odd Y of bit length W reaches an even number within W + 4 iterations of G.
- E implies Q01, and so no eventual period 2 for the centre column.
- The Kit proves only that Q01 implies *some* uniform bound (`uniform_death`, existential). It states no explicit one.

## Questions
- Wesley: prove a stripe-front descent. The Kit's stripe-prefix families show death within 31 steps of arrival (`Rule30ZeroBoundaryStripeFamily`, `Rule30ZeroBoundaryStripeFour`). Do the record holders all come from those families, and does one potential function bound them all?
- Data: does the record grow like c·W with c < 1? Fit the records from 2⁶⁶ onward, which needs W ≥ 67 (2³³ leaves and up).
- Is Conjecture E already known? Check the Kit's Q01Gap and the stripe-front modules.
