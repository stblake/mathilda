---
source: src/special_functions/fresnel.c
references:
  - "DLMF §7.2(iii) — the Fresnel integral C(z) = Int_0^z cos(π t^2/2) dt."
  - "DLMF §7.12 — the asymptotic expansion of the Fresnel integrals."
---
**Algorithm.** `builtin_fresnelc` handles `FresnelC[z] = Int_0^z cos(π t^2/2) dt` (the π/2-normalized / Wolfram convention), entire and odd. It shares one numeric kernel with `FresnelS`: the pair `(C, S)` is computed together and this builtin returns the `C` component. Exact special values first: `FresnelC[0] = 0`, `FresnelC[±Infinity] = ±1/2`, `FresnelC[±I Infinity] = ±I/2`, `ComplexInfinity`/`Indeterminate -> Indeterminate`. A **numeric real** argument (machine or arbitrary precision) uses the convergent Maclaurin series for small/moderate `|x|` (with `~(π/2)|x|^2/ln2` guard bits to absorb the `~e^{(π/2)|z|^2}` partial-sum cancellation), or the **asymptotic expansion** DLMF 7.12 (`C(x) = 1/2 + f(x) sin(π x^2/2) - g(x) cos(π x^2/2)`, summed to optimal truncation) for large `|x|` — the asymptotic constant `1/2` holds only in a sector around the real axis (Stokes), so it is used for real inputs only. A **complex** argument always uses the convergent paired `A/B` series (`A = C + iS`, `B = C - iS`) in the shared `ncpx` toolkit (correct everywhere). A symbolic negative-leading argument folds by oddness. A `USE_MPFR=0` build uses a double-complex `A/B` series.

**Data structures.** `Expr`; the shared complex-MPFR toolkit `ncpx` (`numeric_complex.h`), folding by oddness to `Re z >= 0`. The ND kernel is a real `REG_U` registration (`NDKU_FresnelC`, `ndk_FresnelC_r` → `sf_machine_fresnel_c` in `src/special_functions/sf_machine.c`): element-wise over a packed or visible real `NDArray`. `Compile[]` lowers `FresnelC` at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes.

**Complexity / limits.** `O(1)` per element at machine precision; MPFR term count and guard bits scale with `|z|^2` and precision. Entire function (no branch cuts). `D[FresnelC[x], x] = Cos[π x^2/2]`. Symbolic arguments stay symbolic. Attributes: `Listable`, `NumericFunction`, `Protected`.
