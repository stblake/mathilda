---
source: src/special_functions/beta.c
references:
  - "DLMF §5.12 — the Euler beta function B(a,b) = Gamma(a) Gamma(b)/Gamma(a+b)."
  - "DLMF §8.17 — the incomplete beta function B_z(a,b) and its relation to the Gauss hypergeometric function."
---
**Algorithm.** `builtin_beta` assembles the beta function from primitives that already exist rather than re-deriving the transcendental machinery. `Beta[a, b] = Gamma(a) Gamma(b) / Gamma(a+b)` is reduced to a gamma ratio and handed back to the evaluator (inheriting `Gamma`'s exact integer/half-integer/rational and machine/arbitrary-precision/complex paths). Pole structure on the integer lattice is classified up front: counting `p = [a∈Z≤0] + [b∈Z≤0] - [a+b∈Z≤0]` gives a surviving pole (`ComplexInfinity`), a cancelling pair (a finite `(-1)^B k!(B-1)!/m!` limit), or a denominator-only pole (`0`). When one argument is a positive integer `n` the ratio collapses exactly to `(n-1)!/Pochhammer[b, n]` (so `Beta[3, 1/3] = 27/14`). The incomplete form `Beta[z, a, b] = z^a/a · 2F1(a, 1-b; a+1; z)` routes through `Hypergeometric2F1` (terminating to a closed form when `b` is a positive integer, numeric otherwise); `Beta[0, a, b] = 0`, `Beta[1, a, b] = Beta[a, b]`; and `Beta[z0, z1, a, b] = Beta[z1, a, b] - Beta[z0, a, b]`. A purely symbolic two-argument call stays unevaluated, matching the Wolfram Language.

**Data structures.** `Expr` trees driven through `eval_and_free`. The ND kernel is a binary `REG_B` registration (`NDKB_Beta`): over real arrays it computes `tgamma(a) tgamma(b) / tgamma(a+b)` (libm), declining complex buffers to the `List` path. `Compile[]` lowers `Beta` (e.g. `Beta[2, x]`) at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes via that kernel.

**Complexity / limits.** The two-argument reduction is as fast as the underlying `Gamma`/`Pochhammer`; the incomplete form costs a `2F1` evaluation. Symbolic arguments keep the call inert; non-cancelling gamma poles give `ComplexInfinity`. The four-argument generalized form reduces only when every argument is numeric (symbolic differentiation is handled in `deriv.c`). Attributes: `Listable`, `NumericFunction`, `Protected`.
