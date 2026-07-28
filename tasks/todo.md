# Limit: oscillatory normal form (proving limits do not exist)

(Previous task — Compile[] M5 — archived to `plans/COMPILE_M5_TODO.md`.)

Motivating case:

    Limit[(Cos[x^2]/x^2 - Cos[(x+1)^2]/(x+1)^2)/(1/x^3), x -> Infinity]

Asymptotically `2 x Sin[x^2+x+1/2] Sin[x+1/2] + O(1)` — an unbounded
oscillation with *no* single dominant term, so none of the existing layers
(squeeze envelope, Series, Gruntz, L'Hospital) can touch it.

## The general class

Every existing layer treats an oscillation as an opaque "bounded head". The
generalisation is to put the whole expression into an **oscillatory normal
form** at an infinite limit point:

    f(x) = c_0(x) + SUM_j c_j(x) E^(I theta_j(x))

with the `theta_j` pairwise-distinct real phases carrying **no constant term**
(constants are absorbed into the amplitude) and the amplitudes `c_j`
oscillation-free. `TrigToExp` + `Expand` produces exactly this.

Distinct phases are asymptotically orthogonal, so the normal form is a
*decision* form: nothing can cancel between groups. The verdicts below are
theorems, not heuristics.

## Decision rules

Let `S` = groups with a non-constant phase, `gamma_j = lim |c_j|`.

- **R1 (squeeze).** All `gamma_j = 0` (j in S) ⟹ `lim f = lim c_0`, since
  `|f - c_0| <= SUM |c_j| -> 0`.
- **R3 (dominated oscillation).** Every `|c_j| = o(|c_0|)` and `c_0 -> ±oo`
  ⟹ `lim f = lim c_0`.
- **R0 (strictly dominant oscillation, IVT).** One phase group (with its
  conjugate mate) strictly dominates every other group, its phase `-> ±oo`
  continuously, and `arg c_j` is bounded. Then `f = c_0 + A cos(theta+phi) +
  o(A)` and IVT hands us two sequences with distinct limits ⟹ **no limit**.
  No polynomial restriction on the phase, so this covers `E^x Cos[x]`,
  `x^5 Cos[x]`, `Sin[Log[x]]`.
- **R2 (mean / mean-square, Weyl).** Phases are real polynomials of degree
  `>= 1` with numeric coefficients. Then

      (1/T) INT_0^T |f|^2  ~  SUM_j (1/T) INT_0^T |c_j|^2
      (1/T) INT_0^T f      ~  (1/T) INT_0^T c_0

  (cross terms die by van der Corput: every phase *difference* is a
  non-constant polynomial because constant terms were stripped). If `f -> L`
  finite then `L = lim c_0` and `|L|^2 = |lim c_0|^2 + SUM gamma_j^2`, forcing
  every `gamma_j = 0`. So **some `gamma_j != 0` ⟹ no finite limit**. `±oo` is
  excluded by the window mean `(1/T) INT_T^2T f`, which stays bounded when
  `|c_j| = O(x^deg theta_j)` — or trivially when `f` is bounded.

Verdict for R0/R2 is `Indeterminate`, matching Mathilda's existing convention
(`limit.h`: "Indeterminate -- provably no limit").

## Tasks

- [x] Probe existing behaviour; confirm `TrigToExp`/`Expand`/`PolynomialQ`/
      `Exponent`/`PossibleZeroQ` suffice to build the normal form.
- [x] `src/calculus/limit_osc.{c,h}` — normal form + the four rules.
- [x] Wire into the `compute_limit` cascade (before `layer2_series`) and add
      `Method -> "Oscillatory"`.
- [x] Reduce a finite limit point to `+Infinity` via `x = a ± 1/t`.
- [x] Tests: `tests/test_limit_oscillatory.c`.
- [x] Docs: `docs/spec/builtins/calculus.md` + weekly changelog.
- [x] Full regression run of the limit/calculus suites.

## Review

**What shipped.** `src/calculus/limit_osc.{c,h}` (~840 lines) plus a ~90-line
hookup in `limit.c`. The layer runs after the cheap squeeze envelope and before
`Series` — Series has no expansion at infinity for `Sin[x^2]` and would either
fail or fold an oscillation into a spurious leading term.

**Three implementation traps.**

1. `TrigToExp` returns **un-flattened** `Times[c, Times[x, E^(I x)]]`. The
   evaluator flattens `Times` when it evaluates, but an un-re-evaluated builtin
   result violates the invariant. A factor collector that does not recurse
   through nested same-head nodes silently drops the `E`-factor — `Sin[x]`
   worked, `x Sin[x]` did not. `collect_parts` now recurses.
2. Compare `|c|^2 = c conj(c)`, never `Abs[c]`: `Abs` stays inert on symbolic
   arguments so its limit is undecidable, while `c conj(c)` expands to a real
   rational expression. Conjugation under "all symbols are real" is just
   negating every literal `Complex[a, b]` imaginary part — exact, and unlike
   `Conjugate[]` it actually reduces.
3. Take square roots on the *limits*, not the expressions. R3 needs
   `SUM |c_j| / |c_0| < 1`; computing `lim |c_j|^2/|c_0|^2` per group and
   `Sqrt`-ing those *numbers* decides `x^2 (2 + Cos[x]) -> Infinity` (ratio
   1/2), which the term-by-term `o(c_0)` test I started with could not.

**Three pre-existing wrong answers fell out** (all verified against a stashed
`limit.c`, so none of them are regressions from this work):

- `Limit[E^(I x)/x, x -> Infinity]` gave `E^DirectedInfinity[I]` — the `1/x`
  swallowed. `exp_of_limit` now refuses a folded value that still carries an
  infinity.
- `Limit[Cos[1/x] - Cos[1/x + 1], x -> 0]` gave `0`: substitution turns both
  terms into `Cos[ComplexInfinity]` and the `Plus` cancels them. Both
  substitution fast paths now refuse when an *inner argument* diverges at the
  point — which is exactly discontinuity, where substitution was never a valid
  limit.
- One assertion in `test_limit.c` recorded an abstention that is now a proof
  (`Sin[x^2] + Cos[x]`).

**Verification.** `limit_oscillatory_tests` (9 groups), plus `limit_tests`,
`limit_assumptions_tests`, `gruntz_tests`, `gruntz_stress_tests`,
`series_tests`, `nlimit_tests`, `nseries_tests`, `residue_tests` and the full
`tests/build` sweep. Valgrind on the new path: byte-identical leak totals to
the `Expand[TrigToExp[...]]` baseline (13,376 B / 418 blocks, all dyld/objc).
`make check-c99` clean.

**Known gaps (honest abstentions, documented in the module header).**
`x^2 (1 + Cos[x])` — envelope exactly equal to the smooth part; `x E^(I x)` —
unmated and unbounded, where `ComplexInfinity` may be the intended answer;
`Tan`/`Sec`/`Csc`, whose `TrigToExp` image leaves an exponential in a
denominator; phases that neither diverge nor are polynomial (`Sin[x] +
Sin[x + 1/x]`); symbolic amplitudes (`a Sin[x]`, since `a = 0` has limit 0).

**Worth a separate look.** `TrigToExp` (and therefore `Expand` of its result)
leaving a nested `Times` is a canonicalisation bug independent of `Limit`; it
will bite anything that walks factors structurally.
