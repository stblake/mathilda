---
source: src/special_functions/fresnel.c
references:
  - "DLMF §7.2(iii) — the Fresnel integral S(z) = Int_0^z sin(π t^2/2) dt."
  - "DLMF §7.12 — the asymptotic expansion of the Fresnel integrals."
---
**Algorithm.** `builtin_fresnels` handles `FresnelS[z] = Int_0^z sin(π t^2/2) dt` (the π/2-normalized / Wolfram convention), entire and odd. It shares the numeric kernel with `FresnelC`: the pair `(C, S)` is computed together and this builtin returns the `S` component. Exact special values first: `FresnelS[0] = 0`, `FresnelS[±Infinity] = ±1/2`, `FresnelS[±I Infinity] = ∓I/2`, `ComplexInfinity`/`Indeterminate -> Indeterminate`. A **numeric real** argument (machine or arbitrary precision) uses the convergent Maclaurin series for small/moderate `|x|` (with `~(π/2)|x|^2/ln2` guard bits), or the **asymptotic expansion** DLMF 7.12 (`S(x) = 1/2 - f(x) cos(π x^2/2) - g(x) sin(π x^2/2)`, optimal truncation) for large `|x|` (real inputs only — the `1/2` constant is sector-dependent by the Stokes phenomenon). A **complex** argument always uses the convergent paired `A/B` series (`A = C + iS`, `B = C - iS`; `S = (A - B)/(2i)`) in the shared `ncpx` toolkit (correct everywhere). A symbolic negative-leading argument folds by oddness. A `USE_MPFR=0` build uses a double-complex `A/B` series.

**Data structures.** `Expr`; the shared complex-MPFR toolkit `ncpx` (`numeric_complex.h`), folding by oddness to `Re z >= 0`. The ND kernel is a real `REG_U` registration (`NDKU_FresnelS`, `ndk_FresnelS_r` → `sf_machine_fresnel_s` in `src/special_functions/sf_machine.c`): element-wise over a packed or visible real `NDArray`. `Compile[]` lowers `FresnelS` at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes.

**Complexity / limits.** `O(1)` per element at machine precision; MPFR term count and guard bits scale with `|z|^2` and precision. Entire function (no branch cuts). `D[FresnelS[x], x] = Sin[π x^2/2]`. Symbolic arguments stay symbolic. Attributes: `Listable`, `NumericFunction`, `Protected`.
