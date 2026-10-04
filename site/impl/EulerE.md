---
source: src/special_functions/eulere.c
references:
  - "DLMF §24.2, §24.4 — Euler numbers and polynomials: generating function 2 e^{xt}/(e^t+1), E_n = 2^n E_n(1/2), and the recurrence for E_{2m}."
---
**Algorithm.** `builtin_eulere` serves `EulerE[n]` (the Euler number `E_n`) and `EulerE[n, x]` (the Euler polynomial `E_n(x)`). The numbers come from the recurrence `E_0 = 1`, `E_{2m} = -Sum_{k=0}^{m-1} C(2m,2k) E_{2k}`, with odd indices identically zero, kept as exact GMP integers in a lazily-grown, process-lifetime cache. Layered dispatch: an exact non-negative integer `n` returns the exact integer `E_n`; an inexact integer-valued `n` returns it numericalised (Real/MPFR); `EulerE[n, x]` builds the degree-`n` polynomial in monomial form, `E_n(x) = Sum_i c_i x^i` with `c_i = S_i / 2^{n-i}` and the all-integer inner sum `S_i = Sum_{j=i}^{n} (-1)^{j-i} C(n,j) C(j,i) E_{n-j}` (exact rational coefficients, evaluated once); and a symbolic-`n` special case `EulerE[n, 1/2] -> 2^-n EulerE[n]`. Everything else stays symbolic. Degenerate numeric results (e.g. `EulerE[200.]`, overflowing a `double`) are promoted to an extended-exponent MPFR real.

**Data structures.** `Expr`; GMP `mpz_t` cache for the numbers and `mpq_t` for polynomial coefficients, assembled into canonical `Integer`/`BigInt`/`Rational` leaves. Caps: `EULER_NUMBER_CAP = 5000`, `EULER_POLY_CAP = 1000`. No ND kernel and no `Compile[]` lowering (`CompileDiagnostics` reports `Compiled -> False`): `EulerE` is an integer-indexed exact table, not an element-wise machine-real function. Diagnostics route through `mth_message` (`argt`).

**Complexity / limits.** The number recurrence is `O(n^2)` in growing big-integers; the polynomial is a double sum of exact integer binomials. Beyond the caps the call stays symbolic. Attributes: `Listable`, `Protected`.
