---
source: src/special_functions/cosintegral.c
references:
  - "DLMF §6.2.11, §6.6.5 — the cosine integral Ci(z) and its convergent Maclaurin series."
  - "DLMF §6.12.2 — the asymptotic expansion Ci(z) ~ sin(z) f(z) - cos(z) g(z)."
---
**Algorithm.** `builtin_cosintegral` handles `CosIntegral[z] = Ci(z)`, which has a logarithmic singularity at `0` and a branch cut along `(-Infinity, 0]` (not entire, not odd). Exact special values first: `Ci[0] = -Infinity`, `Ci[Infinity] = 0`, `Ci[-Infinity] = I Pi`, `Ci[±I Infinity] = Infinity`, `ComplexInfinity`/`Indeterminate -> Indeterminate`. A **numeric real or complex** argument routes to the MPFR kernel: for moderate `|z|` the **convergent Maclaurin series** DLMF 6.6.5 (`Ci(z) = γ + Log(z) + Sum_{k>=1} (-1)^k z^{2k}/(2k(2k)!)`, with `|z|/ln2` guard bits), whose principal `Log(z)` supplies the correct `±iπ` jump across the cut with no folding; for large `|z|` the **asymptotic expansion** DLMF 6.12.4 (`sin(z) f(z) - cos(z) g(z)`) plus a piecewise Stokes/reflection constant (`0`, `±π/2`, or `±π`) restoring the principal branch (a negative real `x` gives the from-above value `Complex[Ci(|x|), Pi]`). The complex path uses the shared `ncpx` toolkit. A `USE_MPFR=0` build uses a machine-double series/asymptotic.

**Data structures.** `Expr`; the shared complex-MPFR toolkit `ncpx` (`numeric_complex.h`). The ND kernel is a real `REG_U` registration (`NDKU_CosIntegral`): real buffers via `ndk_CosIntegral_r` → `sf_machine_ci`, complex buffers via `cosintegral_machine_complex` (the file's own machine series, with a shared `sf_series_usable` cancellation gate that declines rather than returns garbage). `Compile[]` lowers `CosIntegral` at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes.

**Complexity / limits.** `O(1)` per element at machine precision; MPFR term count and guard bits scale with `|z|` and precision. Logarithmic singularity at `0`, branch cut on the negative real axis. Symbolic arguments stay symbolic. Attributes: `Listable`, `NumericFunction`, `Protected`.
