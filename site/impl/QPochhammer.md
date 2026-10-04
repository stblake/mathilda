---
source: src/special_functions/qpochhammer.c
references:
  - "DLMF §17.2 — q-series and the q-Pochhammer symbol (q-shifted factorial)."
---
**Algorithm.** `builtin_qpochhammer` evaluates the q-Pochhammer symbol. The
three-argument finite form `QPochhammer[a, q, n] = prod_{k=0}^{n-1} (1 - a q^k)`:
for a non-negative integer `n` the product is built as an `Expr` and handed to
the evaluator, which reduces it exactly for exact `a, q` and at machine / MPFR
precision for inexact `a, q`; `n = 0 -> 1`; a symbolic or non-integer `n` is
left unevaluated (which is exactly what `Product` relies on to emit the closed
form). The two-argument infinite form `QPochhammer[a, q] = (a; q)_inf` is
evaluated only for machine-real `a, q` with `|q| < 1` **and** at least one
inexact operand, by accumulating `(1 - a q^k)` in `double` until the factor
falls below machine epsilon; an all-exact or `|q| >= 1` input stays symbolic.

**Data structures.** `Expr`; the finite form goes through the ordinary evaluator
(GMP for exact, libm/MPFR for inexact); the infinite form is a `double`
accumulation only. ND: binary kernel
`NDK_BIN2(QPochhammer, sf_machine_qpochhammer)`, registered `REG_B`, so
`packed_aware`. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** The finite form is `O(n)` factors; the infinite form is
**machine `double` only** (no MPFR kernel — so `N[QPochhammer[a, q], p]` returns
a machine-precision value regardless of `p`) and requires `|q| < 1`, summing up
to 100000 factors. `Compile[]` lowers at both scalar and rank-1 array shapes
(`Compiled -> True`).
