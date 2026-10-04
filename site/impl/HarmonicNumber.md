---
source: src/special_functions/harmonicnumber.c
references:
  - "DLMF §25.11 — the Hurwitz zeta function and the relation H_n^{(r)} = ζ(r) - ζ(r, n+1)."
  - "DLMF §5.15 — the digamma function (H_n = γ + ψ(n+1))."
---
**Algorithm.** `builtin_harmonicnumber` serves `HarmonicNumber[n] = Sum_{i=1}^n 1/i` and the generalized `HarmonicNumber[n, r] = Sum_{i=1}^n 1/i^r`. Rather than carry bespoke numeric kernels it reduces to existing primitives and lets the evaluator finish: an exact **non-negative integer `n`** (within `HN_EXPAND_CAP = 100000`) expands to the explicit finite sum `Sum_{i=1}^n i^{-r}` (an exact rational for integer `r`, an explicit `Plus` for symbolic/complex `r`); **`n -> Infinity`** gives `Zeta[r]`; a **non-positive integer `r = -m`** gives the Faulhaber polynomial in `n` built from `BernoulliB`; and an **inexact / numericizable** argument uses the analytic identity `H_n^{(r)} = Zeta[r] - Zeta[r, n+1]` (or, for `r == 1`, `EulerGamma + PolyGamma[0, n+1]`) wrapped in `N[…]` at the input precision — so arbitrary precision and complex arguments pass straight through `Zeta`/`PolyGamma`. A `numericizable` guard keeps `HarmonicNumber[x, 2.5]` symbolic in a free symbol `x`. Everything else stays symbolic.

**Data structures.** `Expr` trees driven through `eval_and_free`; GMP for the integer/exact paths. The ND kernel is a real `REG_U` registration (`NDKU_HarmonicNumber`, `ndk_HarmonicNumber_r` → `sf_machine_harmonic`, which evaluates `γ + ψ(x+1)` via the machine digamma): element-wise over a packed or visible real `NDArray`. `Compile[]` lowers `HarmonicNumber` at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes.

**Complexity / limits.** The finite-sum expansion is `O(n)` and capped at `n <= 100000` (an exact integer argument is never numerically contaminated, so beyond the cap it simply stays symbolic — there is no cheaper fallback). Numeric reduction cost is that of `Zeta`/`PolyGamma`. Diagnostics route through `mth_message` (`argt`). Attributes: `Listable`, `NumericFunction`, `Protected`.
