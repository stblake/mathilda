---
source: src/special_functions/logintegral.c
references:
  - "DLMF §6.2 — the logarithmic integral li."
---
**Algorithm.** `builtin_logintegral` evaluates `li(z)` through the identity
`li(z) = Ei(Log z)`, reusing `ExpIntegralEi`'s numeric stack rather than
duplicating it; the principal `Log` supplies the `±i Pi` jump that places the
branch cut on `(-Infinity, +1)`. Exact special values: `0 -> 0`,
`1 -> -Infinity`, `Infinity -> Infinity`, `ComplexInfinity`/`Indeterminate ->
Indeterminate`. Exact non-special numbers (`li[2]`, `li[1/2]`) stay symbolic,
matching the Wolfram Language; only an inexact argument (or explicit `N[...]`)
evaluates, by building `ExpIntegralEi[Log[z]]` and evaluating it — a defensive
check drops the result and stays symbolic if the composition did not reduce. The
machine kernel `logintegral_machine_complex` is `Ei(clog z)`, declining at
`z = 0` (Indeterminate) and `z = 1` (pole).

**Data structures.** `Expr`; reuses the `ExpIntegralEi` kernels (`mpfr_eint` /
the real and complex convergent series); `double complex` machine kernel. ND:
unary kernel `NDKU_LogIntegral = { logintegral_machine_complex,
ndk_LogIntegral_r, ... }` (the real kernel is `sf_machine_li`), registered
`REG_U`, so `packed_aware`. Attributes: `Listable`, `NumericFunction`,
`Protected`.

**Complexity / limits.** Cost is one `ExpIntegralEi` evaluation. The principal
`Log` puts the cut on the negative real axis, so `li[-1.]` is complex
(`Ei(i Pi)`). Exact rationals stay symbolic. `Compile[]` lowers at both scalar
and rank-1 array shapes (`Compiled -> True`).
