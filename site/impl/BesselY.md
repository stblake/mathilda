---
source: src/special_functions/bessel.c
references:
  - "DLMF §10.8.1 — the logarithmic series for integer-order Y_n(z) near the origin."
  - "DLMF §10.17.4 — the large-argument asymptotic expansion of Y_nu(z)."
---
**Algorithm.** `builtin_bessely` (via `bessely_two_arg`) handles `BesselY[n, z]`, the second-kind solution, singular at the origin. Exact `z == 0` with classified order gives `Y_0(0) = -Infinity`, `0` (a negative half-odd-integer order, where `cos(nu π) = 0` cancels the divergent term so `Y_{-1/2} = J_{1/2}`), `Indeterminate` (`Re nu = 0`, `nu != 0`), or `ComplexInfinity` otherwise. For a numeric call: integer order and real `z > 0` take the MPFR-native **`mpfr_yn`** fast path (correctly rounded); otherwise the unified core `by_core` routes among (i) the **connection** `Y_nu(z) = (J_nu(z) cos(nu π) - J_{-nu}(z)) / sin(nu π)` for small `|z|`, non-integer order, (ii) the **logarithmic series** DLMF 10.8.1 for small `|z|`, integer order, and (iii) the **asymptotic series** DLMF 10.17.4 (`Y_nu(z) ~ sqrt(2/(πz))[sin(w) A + cos(w) B]`) for large `|z|`, summed to optimal truncation. Half-integer→elementary rewrites live in `src/internal/bessel.m`; `Series`/`D` in `calculus/`. Everything else stays symbolic.

**Data structures.** `Expr`; the shared complex-MPFR toolkit `ncpx` at explicit working precision. The ND kernel is a binary `REG_B` registration (`NDKB_BesselY`): over real arrays with **integer** order it calls libc `yn` (declines non-integer order or any complex operand to the `List` path). `Compile[]` lowers `BesselY[n, z]` at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes via that kernel.

**Complexity / limits.** `mpfr_yn` is `O(1)` at machine precision for integer order; the core's term counts and guard bits scale with `|z|` and precision. Branch cut along the negative real `z` axis. Symbolic arguments stay symbolic. Attributes: `Listable`, `NumericFunction`, `Protected`.
