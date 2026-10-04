---
source: src/special_functions/erfi.c
references:
  - "DLMF §7.2, §7.6.2 — the imaginary error function erfi(z) = erf(iz)/i = (2/sqrt(π)) Int_0^z e^{t^2} dt."
---
**Algorithm.** `builtin_erfi` handles `Erfi[z] = -i erf(iz)`, entire and odd. There is no libm `erfi` and no `mpfr_erfi`, so both numeric kernels are hand-rolled. Exact special values first: `Erfi[0] = 0`, `Erfi[±Infinity] = ±Infinity`, `Erfi[ComplexInfinity] = ComplexInfinity`, and the finite imaginary-axis limits `Erfi[±I Infinity] = ±I`, `Erfi[Indeterminate] = Indeterminate`. A **real** argument (machine or arbitrary precision) uses the all-positive Maclaurin series `erfi(x) = (2/sqrt(π)) Sum x^{2n+1}/(n!(2n+1))` in MPFR — every term shares `x`'s sign, so the partial sums climb monotonically with no cancellation and only a flat 64-bit guard. A **complex** argument (any precision) uses `erfi(z) = -i erf(iz)`, reusing the cancellation-aware `erf` Maclaurin series DLMF 7.6.2 in the file-local `ecx` toolkit with `|z|^2/ln2` guard bits. A symbolic negative-leading argument folds by oddness (`Erfi[-x] = -Erfi[x]`). A `USE_MPFR=0` build falls back to double-precision series.

**Data structures.** `Expr`; file-local `ecx` (pairs of `mpfr_t`) for the complex path, plain MPFR scalars for the real series. The ND kernel is a real `REG_U` registration (`NDKU_Erfi`, `ndk_Erfi_r` → `sf_machine_erfi` in `src/special_functions/sf_machine.c`): element-wise over a packed or visible real `NDArray`. `Compile[]` lowers `Erfi` at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes. Degenerate machine results are promoted via `numeric_promote_result_if_degenerate`.

**Complexity / limits.** `O(1)` per element at machine precision; the real series needs no cancellation guard, the complex series' term count and guard bits scale with `|z|^2` and precision. Entire function (no branch cuts). Symbolic arguments stay symbolic. Attributes: `Listable`, `NumericFunction`, `Protected`.
