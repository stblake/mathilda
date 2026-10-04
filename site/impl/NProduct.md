---
source: src/numerical_calculus/nprod.c
references:
  - "J. B. Keiper, *Numerical computation of infinite products*, Wolfram Research tech. report (1992) — the Exp[NSum[Log]] reduction."
---
**Algorithm.** `builtin_nproduct` (HoldAll) parses the iterator `{i, imin, imax,
di}` and evaluates the product as **`Exp[NSum[Log[f(i)], {i, imin, imax}]]`**
(Keiper 1992): it forms `logbody = Log[f]`, builds an inner `NSum` with the
options mapped across (`NProductFactors -> NSumTerms`, so the default leading
factor count is NSum's 15; `NProductExtraFactors -> NSumExtraTerms`; `Method`,
`WynnDegree`, `VerifyConvergence`, the goals and `WorkingPrecision` passed
through), evaluates it, and returns `Exp[result]`. So the whole NSum engine —
Euler–Maclaurin for monotone factors, Wynn's epsilon otherwise, and its
divergence test — carries over directly; see the `NSum` notes for the summation
machinery. Multidimensional products nest an inner `NProduct` as the factor body
(dependent inner bounds see the outer index via HoldAll). A divergent log-sum
(the inner `NSum` returning `ComplexInfinity`) makes `NProduct` return
`ComplexInfinity`; a non-numeric inner result leaves the product unevaluated.

**Data structures.** No numerical buffers of its own — everything is the inner
`NSum`'s (machine `double _Complex` or split `mpfr_t` real/imag, with the Log
summand drawn through NSum's cached compiled programs). The only `NProduct`-
specific state is the option-remapping and the final `Exp`.

**Complexity / limits.** The cost is one `NSum` of `Log[f]` plus one `Exp`. The
one numerical subtlety is precision: `Exp` amplifies the exponent's *absolute*
error into the product's *relative* error, so the inner sum runs at
`WorkingPrecision + 10` guard digits and the final `Exp` is renormalised back to
the requested precision with `N[·, wdigits]`. Options: `Method` (`Automatic` |
`EulerMaclaurin` | `WynnEpsilon`), `WorkingPrecision`, `NProductFactors`,
`NProductExtraFactors`, `WynnDegree`, `VerifyConvergence` (default `True`; a
divergent product gives `ComplexInfinity`), `AccuracyGoal`, `PrecisionGoal`.
`NProduct` is `HoldAll` and `Protected`.
