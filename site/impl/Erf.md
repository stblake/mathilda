---
source: src/special_functions/erf.c
references:
  - "DLMF §7.2, §7.6.2 — the error function and its Maclaurin series erf(z) = (2/sqrt(π)) e^{-z^2} Sum_n t_n."
---
**Algorithm.** `builtin_erf` handles `Erf[z]` (and `Erf[z0, z1] = erf(z1) - erf(z0)`). `erf` is entire and odd. Exact special values first: `Erf[0] = 0`, `Erf[±Infinity] = ±1`, `Erf[ComplexInfinity] = ComplexInfinity`, `Erf[±I Infinity] = DirectedInfinity[±I]`, `Erf[Indeterminate] = Indeterminate`. Then by argument kind: a **machine real** uses libm `erf`; an **arbitrary-precision real** uses `mpfr_erf`; a **complex** argument (any precision) uses the cancellation-aware Maclaurin series DLMF 7.6.2 (`erf(z) = (2/sqrt(π)) e^{-z^2} Sum t_n`, `t_0 = z`, `t_n = t_{n-1}(2z^2)/(2n+1)`) evaluated in the file-local complex-MPFR toolkit `ecx` with `|z|^2/ln2` guard bits to absorb the `~e^{|z|^2}` partial-sum cancellation exactly — so even machine-precision complex results carry full accuracy. A symbolic negative-leading-coefficient argument folds by oddness (`Erf[-x] = -Erf[x]`). `Erf[z0, z1]` commits only if both erf calls become concrete. A `USE_MPFR=0` build falls back to a double-complex series. An `Interval` argument routes to `interval_apply_function`.

**Data structures.** `Expr`; file-local `ecx` (pairs of `mpfr_t`, alias-safe, explicit precision) for the complex series. The ND kernel is a real `REG_U` registration (`NDKU_Erf`, `ndk_Erf_r` → libm `erf`): element-wise over a packed or visible real `NDArray`, no complex ND kernel. `Compile[]` lowers `Erf` at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes.

**Complexity / limits.** `O(1)` per element on the real paths; the complex series term count and guard bits scale with `|z|^2` and precision. Entire function (no branch cuts). Symbolic arguments stay symbolic. Diagnostics route through `mth_message`. Attributes: `Listable`, `NumericFunction`, `Protected`.
