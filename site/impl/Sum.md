---
references:
  - "R. W. Gosper, \"Decision procedure for indefinite hypergeometric summation\", Proc. Natl. Acad. Sci. USA 75 (1978)."
  - "M. Petkovšek, H. S. Wilf, D. Zeilberger, *A = B* (A K Peters, 1996)."
source: src/sum/sum_gosper.c
---
**Algorithm.** `Sum` is `HoldAll` so the iterator variable and bounds are not
prematurely evaluated. The dispatcher `builtin_sum` (src/sum/sum.c) strips
trailing options, rewrites multiple iterators `Sum[f, s1, ..., sk]` into nested
single-spec sums (outer bounds may depend on inner variables), and then:

- **Finite explicit expansion** — when a range resolves to a finite integer
  span, or the spec iterates an explicit list, it binds the variable and folds
  the evaluated terms with `Plus` (runaway-guarded by `SUM_MAX_FINITE_TERMS`).
- **Closed-form cascade** — for symbolic bounds, `Infinity`, or the indefinite
  form `Sum[f, i]`, it runs the context-qualified sub-algorithms in order:
  `Sum\`Polynomial` (Faulhaber-style polynomial/telescoping sums,
  src/sum/sum_polynomial.c), `Sum\`Geometric` (geometric/exponential terms,
  src/sum/sum_geometric.c), then `Sum\`Gosper` (src/sum/sum_gosper.c). Each
  stage returns the closed form or comes back unevaluated to fall through; if
  all fall through, `Sum[...]` is returned held. A `Method -> "..."` option
  selects a single strict stage.

The key stage is **Gosper's algorithm** for indefinite hypergeometric
summation. Given a term `t(i)` whose ratio `r(i) = t(i+1)/t(i)` is rational
(`Simplify` reduces factorial ratios, `Together` yields polynomials `a, b`),
it computes the **Gosper–Petkovšek normal form** `r = (a/b)·(c(i+1)/c(i))` with
`gcd(a(i), b(i+h)) = 1` for all integers `h ≥ 0`, found via the dispersion set
(`GOSPER_DISPERSION_MAX`-capped) and GCD peeling. It then solves
`a(i)·x(i+1) − b(i−1)·x(i) = c(i)` for a polynomial `x` by undetermined
coefficients (`SolveAlways`); no solution proves `t` is not Gosper-summable.
The antidifference is `F(i) = (b(i−1)/c(i))·x(i)·t(i)`; the definite finite sum
is `F(imax+1) − F(imin)`. The output stays elementary (`R(i)·t(i)`).

**Data structures.** `Expr*` trees throughout; the Gosper stage builds on
Mathilda's polynomial builtins (`Expand`, `PolynomialGCD`, `PolynomialQuotient`,
degree queries) and `SolveAlways`.

**Complexity / limits.** Finite expansion is linear in the number of terms.
Gosper's procedure is a complete decision procedure for hypergeometric
antidifferences, but only for hypergeometric terms; non-hypergeometric or
non-summable inputs fall through to the held form. No creative-telescoping
(Zeilberger) for parametric/definite hypergeometric sums.
