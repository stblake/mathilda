# Hypergeometric0F1

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Hypergeometric0F1[b, z]`**

is the confluent hypergeometric 0F1, equal to HypergeometricPFQ\[{}, {b}, z\].

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Hypergeometric0F1[1/2, z]
Out[1]= Cosh[2 Sqrt[z]]
```

### Applications (3)

```mathematica
In[2]:= Hypergeometric0F1[1/2, z]
Out[2]= Cosh[2 Sqrt[z]]

In[3]:= Hypergeometric0F1[3/2, z]
Out[3]= (1/2 Sinh[2 Sqrt[z]])/Sqrt[z]

In[4]:= N[Hypergeometric0F1[3/2, -(Pi^2/16)]*Pi/2, 40]
Out[4]= 1.0 + 2.546305796115700172884063021572241471189e-60*I
```

## Implementation notes

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

**Attributes:** `NumericFunction`, `Protected`.

## References

- DLMF §16.2 — generalized hypergeometric series; 0F1 is the confluent limit.
- Source: [`src/special_functions/hypergeopfq.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/hypergeopfq.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_hypergeopfq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergeopfq.c)

## Notes & additional examples

### Notes

`Hypergeometric0F1[b, z]` is the confluent limit `HypergeometricPFQ[{}, {b}, z]`. It converges for all `z` and underlies the Bessel functions: `0F1[1/2, z] = Cosh[2 Sqrt[z]]` and `0F1[3/2, z] = Sinh[2 Sqrt[z]]/(2 Sqrt[z])`. The tiny imaginary residue in the last result is numerical noise from the radical of a negative argument.
