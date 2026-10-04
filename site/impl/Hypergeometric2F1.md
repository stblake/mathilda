---
source: src/special_functions/hypergeopfq.c
references:
  - "DLMF §15 — the Gauss hypergeometric function 2F1."
---
**Algorithm.** `builtin_hypergeometric_2f1` is a convenience head:
`Hypergeometric2F1[a, b, c, z]` rewrites to `HypergeometricPFQ[{a, b}, {c}, z]`
and inherits that engine — generic cancellation, `z == 0 -> 1`, termination to a
polynomial when `a` or `b` is a non-positive integer, and the elementary
reductions kept in `try_reduce` (`2F1(1,1;2;z) = -Log[1-z]/z`; the
central-binomial `ArcSin` forms `2F1(1,1;1/2;z)`, `2F1(1,1;3/2;z)`,
`2F1(2,1;3/2;z)`). Numeric is the direct series, valid for `|z| < 1`
(`p = q + 1`); `|z| >= 1` stays unevaluated (no analytic continuation of the
Gauss series here). `LegendreP`/`LegendreQ` numeric paths build on this head.

**Data structures.** `Expr`; the `PFQ` machine `double complex` / MPFR `cpx_t`
summation; `sf_machine_hyper2f1` (wrapping `sf_machine_pfq`) is the shared
machine kernel. ND: N-ary kernel `NDKN_Hypergeometric2F1` (arity 4), element-wise
over the `z` buffer; on `AWARE` with the family. Attributes: `NumericFunction`,
`Protected` (not `Listable`).

**Complexity / limits.** Numeric only for `|z| < 1`; term caps as in
`HypergeometricPFQ`. `Compile[]` lowers at scalar shape (`Compiled -> True`) but
**not** at rank-1 array shape (`Compiled -> False`); the ND kernel still serves
the packed / visible-`NDArray` fast path at the REPL.
