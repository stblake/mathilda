---
source: src/special_functions/bernoullib.c
references:
  - "DLMF §24.2, §24.4 — Bernoulli numbers and polynomials: generating function, the recurrence, and B_n(x) = Sum_j C(n,j) B_{n-j} x^j."
---
**Algorithm.** `builtin_bernoullib` serves `BernoulliB[n]` (the Bernoulli number `B_n`) and `BernoulliB[n, x]` (the Bernoulli polynomial `B_n(x)`). The numbers are computed by the recurrence `B_0 = 1`, `B_m = -1/(m+1) Sum_{k=0}^{m-1} C(m+1,k) B_k`, kept as exact GMP rationals in a lazily-grown, process-lifetime cache; odd indices above 1 are returned as an exact `0`. Argument handling is layered: an exact non-negative integer `n` returns the exact rational `B_n`; an inexact integer-valued `n` (Real/MPFR) returns that rational numericalised to the input precision; `BernoulliB[n, x]` builds `B_n(x) = Sum_{j=0}^{n} C(n,j) B_{n-j} x^j` with exact rational coefficients and evaluates it once (exact `x` stays exact, inexact `x` evaluates numerically, zero coefficients skipped). Everything else — negative, non-integer, or out-of-range `n` — stays symbolic. Degenerate numeric results (e.g. `BernoulliB[300.]`, which overflows a `double`) are promoted to an extended-exponent MPFR real.

**Data structures.** `Expr`; GMP `mpq_t` cache for the numbers, built into canonical `Integer`/`BigInt`/`Rational` leaves. Caps guard runaway: `BERNOULLI_NUMBER_CAP = 5000`, `BERNOULLI_POLY_CAP = 1000`. No ND kernel and no `Compile[]` lowering (`CompileDiagnostics` reports `Compiled -> False`): `BernoulliB` indexes an integer-indexed exact-rational table, not an element-wise machine-real function. Diagnostics route through `mth_message` (`argt`).

**Complexity / limits.** The number recurrence is `O(n^2)` in growing big-rationals; the polynomial expansion is `O(n)` terms with exact rational coefficients. Beyond the caps the call stays symbolic. Attributes: `Listable`, `Protected`.
