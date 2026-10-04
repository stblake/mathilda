---
source: src/special_functions/hypergeopfq.c
references:
  - "DLMF §13 — confluent hypergeometric functions (Kummer's 1F1)."
---
**Algorithm.** `builtin_hypergeometric_1f1` is a convenience head:
`Hypergeometric1F1[a, b, z]` rewrites to `HypergeometricPFQ[{a}, {b}, z]` and
inherits that engine — generic cancellation (`a == b` collapses to `Exp[z]`),
`z == 0 -> 1`, termination to a polynomial when `a` is a non-positive integer,
the reduction `1F1(1;2;z) = (E^z - 1)/z`, and the direct series sum. Numeric is
a machine `double complex` sum or, above machine precision, an MPFR `(re,im)`
sum; because the confluent series alternates for negative real `z` (peak term
`~ e^|z|`), the lost-bits re-sum in MPFR keeps a machine-precision answer
accurate. Being `1F1` (`p = 1 = q`), the series is entire, so a finite `z`
always yields a value in the convergent regime; a machine overflow (e.g.
`Hypergeometric1F1[1., 2., 900.]`) promotes the whole call to MPFR.

**Data structures.** `Expr`; the `PFQ` machine `double complex` / MPFR `cpx_t`
summation; `sf_machine_hyper1f1` (wrapping `sf_machine_pfq`) is the shared
machine kernel. ND: N-ary kernel `NDKN_Hypergeometric1F1` (arity 3), element-wise
over the `z` buffer; on `AWARE` with the family. Attributes: `NumericFunction`,
`Protected` (not `Listable`).

**Complexity / limits.** Entire in `z`; term caps as in `HypergeometricPFQ`.
`Compile[]` lowers at scalar shape (`Compiled -> True`) but **not** at rank-1
array shape (`Compiled -> False`); the ND kernel still serves the packed /
visible-`NDArray` fast path at the REPL.
