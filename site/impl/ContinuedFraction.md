---
source: src/contfrac.c
---
**Algorithm.** `builtin_continued_fraction` computes the *simple* continued-fraction expansion, dispatching on input regime: (1) **exact rationals** use the Euclidean algorithm with floor division (`cf_rational`), producing the canonical terminating form (last term `>= 2`); (2) **`Sqrt[D]`** for a non-square positive integer `D` uses the classic periodic-surd recurrence `m'=a q-m, q'=(D-m'^2)/q, a'=floor((a0+m')/q')`, detecting the period when the `(m, q)` state first repeats (`cf_sqrt_period`), and returns `{a0, {period...}}` (or unrolls to `n` terms); (3) **inexact reals** (machine or MPFR) extract terms by repeated reciprocation while tracking absolute uncertainty (`|x| 2^-prec` plus 64 guard bits), stopping when the integer part is no longer determined by available precision (`cf_inexact`); (4) **exact symbolic reals with explicit `n`** (Pi, etc.) are numericised via `N[x, digits]` at doubling precision until `n` terms stabilise (`cf_exact_numeric`).

**Data structures.** Terms accumulate in a growable `TermVec` of GMP `mpz_t`; the MPFR path uses `mpfr_t` working registers; output is a `List` of integers (with a nested `List` for the repeating block in the unbounded surd case).

**Complexity / limits.** Rationals terminate in `O(log)` Euclidean steps. General quadratic irrationals beyond bare `Sqrt[D]` (e.g. `(1+Sqrt[5])/2`) are not recognised symbolically — supply an explicit `n` to use the numeric path. Inexact inputs stop at the precision floor.
