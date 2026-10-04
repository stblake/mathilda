---
source: src/special_functions/hypergeopfq.c
references:
  - "DLMF §16.2 — generalized hypergeometric series; 0F1 is the confluent limit."
---
**Algorithm.** `builtin_hypergeometric_0f1` is a convenience head:
`Hypergeometric0F1[b, z]` rewrites to `HypergeometricPFQ[{}, {b}, z]` (via
`rebuild_eval`) and inherits that engine entirely — generic cancellation,
`z == 0 -> 1`, termination at a non-positive-integer upper parameter (none here,
`p = 0`), the elementary reductions `0F1(;1/2;z) = Cosh[2 Sqrt z]` and
`0F1(;3/2;z) = Sinh[2 Sqrt z]/(2 Sqrt z)`, and the direct series sum
(machine `double complex`, or MPFR `(re,im)` when precision > 53). Since `p = 0`,
`q = 1` the series is entire, so a numeric value is produced for all finite `z`.

**Data structures.** `Expr`; the `PFQ` machine `double complex` / MPFR `cpx_t`
summation. ND: binary kernel `NDK_BIN2(Hypergeometric0F1, sf_machine_hyper0f1)`
(`sf_machine_hyper0f1(a, z)` wraps `sf_machine_pfq`), element-wise over the `z`
buffer; on `AWARE` with the rest of the family. Attributes: `NumericFunction`,
`Protected` (not `Listable`).

**Complexity / limits.** Entire in `z`; cost is the series term count (capped as
in `HypergeometricPFQ`). `Compile[]` lowers at both scalar and rank-1 array
shapes (`Compiled -> True`), because it is a *binary* ND kernel — unlike the
`1F1`/`2F1`/`PFQ` heads, which lower only at scalar shape.
