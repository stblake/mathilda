---
source: src/poly/algebraicnumberdenominator.c
references:
  - "H. Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138 (Springer, 1993), §4.1–4.2 (algebraic integers; denominators of algebraic numbers)."
  - "The FLINT library (https://flintlib.org), `qqbar` module — minimal polynomial and `fmpz_factor`."
---
**Algorithm.** `builtin_algebraicnumberdenominator` checks arity 1 and delegates
to `flint_qqbar_algebraic_number_denominator`, which computes exactly the
smallest positive integer `d` such that `d·x` is an algebraic integer. It is
**not** `qqbar_denominator` (the leading coefficient `a_n` of the primitive
integer minimal polynomial), which only over-estimates. Writing the monic
minimal polynomial of `d·x`, its `x^i` coefficient is `(a_i/a_n)·d^{n−i}`, so
`d·x` is an algebraic integer iff for every `i < n` the denominator `q_i` of
`a_i/a_n` divides `d^{n−i}`. Since `q_i | a_n`, the minimal `d` divides `a_n`.
The engine factors `a_n` once (`fmpz_factor`) and, per prime `P | a_n`, takes

```
v_P(d) = max_{i<n} ceil( v_P(q_i) / (n−i) ),   v_P(q_i) = max(0, v_P(a_n) − v_P(a_i))
```

(a zero coefficient `a_i` contributes no constraint). Worked example from the
source: `1/5 + Sqrt[2]` has `p = 25x² − 10x − 49`, `a_n = 25`, yet `d = 5`.

A non-constant-algebraic argument gives the engine result `0`, which routes a
`AlgebraicNumberDenominator::nalg` message through `mth_message` (Quiet/Check
funnel) and leaves the expression unevaluated; `−1` (FLINT off) is a silent
decline.

**Data structures.** FLINT `qqbar_t` (minimal polynomial as `fmpz_poly`), an
`fmpz_factor_t` of the leading coefficient, and an `fmpz` valuation scan over the
lower coefficients. The result `d` is returned as an `Integer`/`BigInt` `Expr`.

**Complexity / limits.** One factorisation of the leading coefficient `a_n`
(small in practice) plus an `O(n)` valuation scan. `Listable`, `Protected`;
bounded by the degree cap `QQBAR_DEGREE_CAP = 120`. For any algebraic integer the
answer is `1`.
