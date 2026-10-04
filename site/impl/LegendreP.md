---
source: src/special_functions/legendre.c
references:
  - "DLMF §14 — Legendre and related functions (§14.3 series, §14.7 integer degree)."
---
**Algorithm.** `builtin_legendre_p` handles `LegendreP[n, x]`,
`LegendreP[n, m, x]` and `LegendreP[n, m, a, x]`. For `P_n(x)` an exact integer
order builds the explicit degree-`|n'|` polynomial from the three-term
recurrence `k P_k = (2k-1) x P_{k-1} - (k-1) P_{k-2}` with exact `mpq`
coefficients (using `P_{-1-n} = P_n`, cap `LEG_POLY_CAP = 2000`); `x == 1 -> 1`
for any order. A non-integer order with an inexact argument is evaluated by the
Gauss series `P_n(x) = 2F1(-n, n+1; 1; (1-x)/2)` summed in the `ncpx` MPFR-complex
toolkit (real/complex, machine/arbitrary precision; requires `|(1-x)/2| < 1`).
The associated forms (integer `n`, integer `m >= 0`): type 1 is the Rodrigues
derivative `(-1)^m (1-x^2)^(m/2) d^m/dx^m P_n(x)` (0 when `m > |n'|`); types 2
and 3 are the regularized Gauss polynomial `2F1Reg(-n, n+1, 1-m, (1-x)/2)` times
a `(1±x)^(±m/2)` prefactor. Non-integer/negative `m`, and `|(1-x)/2| >= 1`, stay
symbolic.

**Data structures.** `Expr`; GMP `mpq_t` coefficient arrays for the polynomial;
`ncpx` (`mpfr_t` re/im) for the numeric series. ND: binary kernel
`NDKB_LegendreP = { ndk_LegendreP_c, ... }` (complex, threads `n` and `x`
element-wise), registered `REG_B`, so `packed_aware`. Attributes: `Listable`,
`NumericFunction`, `Protected`.

**Complexity / limits.** The recurrence is `O(n^2)` in growing big rationals,
capped at `n = 2000`; the numeric series costs `O(wp)` terms inside the
convergence disk. Deferred (left symbolic): symbolic `Series`/`D` rules, the
non-integer associated forms, and continuation for `|(1-x)/2| >= 1`. `Compile[]`
lowers at both scalar and rank-1 array shapes (`Compiled -> True`).
