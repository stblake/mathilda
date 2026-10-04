---
source: src/special_functions/erfc.c
references:
  - "DLMF §7.2, §7.6.2 — the complementary error function erfc(z) = 1 - erf(z)."
---
**Algorithm.** `builtin_erfc` handles `Erfc[z] = 1 - erf(z)`, entire with no branch cuts and (unlike `erf`) no symmetry that simplifies `Erfc[-x]`. Exact special values first: `Erfc[0] = 1`, `Erfc[Infinity] = 0`, `Erfc[-Infinity] = 2`, `Erfc[ComplexInfinity] = ComplexInfinity`, `Erfc[±I Infinity] = DirectedInfinity[∓I]` (negated relative to `erf`), `Erfc[Indeterminate] = Indeterminate`. Then by argument kind: a **machine real** uses libm `erfc`; an **arbitrary-precision real** uses `mpfr_erfc` (cancellation-free even for large positive `z`, where `1 - erf(z)` would lose all significance); a **complex** argument (any precision) computes `erf(z)` from the cancellation-aware Maclaurin series DLMF 7.6.2 in the file-local `ecx` toolkit with `|z|^2/ln2` guard bits, then forms the complement `1 - erf(z)` at working precision before rounding (so the guard bits also absorb the `1 - erf` cancellation that dominates for large positive `Re z`). A `USE_MPFR=0` build falls back to a double-complex series. An `Interval` argument routes to `interval_apply_function`.

**Data structures.** `Expr`; file-local `ecx` (pairs of `mpfr_t`) for the complex series. The ND kernel is a real `REG_U` registration (`NDKU_Erfc`, `ndk_Erfc_r` → libm `erfc`): element-wise over a packed or visible real `NDArray`, no complex ND kernel. `Compile[]` lowers `Erfc` at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes. Degenerate machine results are promoted via `numeric_promote_result_if_degenerate`.

**Complexity / limits.** `O(1)` per element on the real paths; the complex series term count and guard bits scale with `|z|^2` and precision. No odd-symmetry fold (`erfc(-x) = 2 - erfc(x)` is left unexpanded). Symbolic arguments stay symbolic. Attributes: `Listable`, `NumericFunction`, `Protected`.
