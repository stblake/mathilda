---
source: src/special_functions/legendre.c
references:
  - "DLMF §14 — Legendre functions of the second kind (§14.3, §14.7)."
---
**Algorithm.** `builtin_legendre_q` handles `LegendreQ[n, x]`,
`LegendreQ[n, m, x]` and `LegendreQ[n, m, a, x]`. Unlike `P_n`, `Q_n` carries a
logarithm even at integer order: for integer `n >= 0` (cap `LEG_POLY_CAP = 2000`,
`n < 0` singular → symbolic) it builds `Q_n(x) = P_n(x) L(x) + v_n(x)` with
`L(x) = (1/2)(Log[1+x] - Log[1-x])` and the polynomial `v_n` from the same
three-term recurrence as `P_n` but seeded `v_0 = 0`, `v_1 = -1`. A non-integer
order emits the exact special value `Q_v(0) = -(Sqrt[Pi]/2) Sin[v Pi/2]
Gamma[(v+1)/2] / Gamma[v/2+1]` for an exact zero argument (which also lets the
origin `Series` fall out of Taylor-via-`D`), and for an inexact argument on the
cut `|x| < 1` evaluates the two Frobenius `2F1` series
`Q_v(0) 2F1(-v/2, (v+1)/2; 1/2; x^2) + Q_v'(0) x 2F1((1-v)/2, (v+2)/2; 3/2; x^2)`
built as an `Expr` over `Hypergeometric2F1`/`Gamma`/`Sin`/`Cos`, so it inherits
their machine, MPFR and complex numerics. Associated `Q_n^m` (all three types,
integer `n, m >= 0`) differentiate `Q_n` in a fresh dummy variable then
substitute.

**Data structures.** `Expr`; GMP `mpq_t` coefficient arrays; the numeric path is
an `Expr` over `2F1`/`Gamma`. ND: binary kernel
`NDKB_LegendreQ = { ndk_LegendreQ_c, ... }`, registered `REG_B`, so
`packed_aware`. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Recurrence `O(n^2)` big rationals, cap `n = 2000`; the
numeric series converges only inside `|x| < 1` (the result is accepted only if
the series actually collapsed to a number, else symbolic). `Compile[]` lowers at
both scalar and rank-1 array shapes (`Compiled -> True`).
