# Hypergeometric1F1

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Hypergeometric1F1[a, b, z]`**

is Kummer's confluent hypergeometric 1F1, equal to HypergeometricPFQ\[{a}, {b}, z\].

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Hypergeometric1F1[a, b, 0]
Out[1]= 1
```

### Applications (4)

```mathematica
In[2]:= Hypergeometric1F1[1, 2, z]
Out[2]= (-1 + E^z)/z

In[3]:= Hypergeometric1F1[-2, 1, z]
Out[3]= 1 - 2 z + 1/2 z^2

In[4]:= N[Hypergeometric1F1[1, 2, 3], 40]
Out[4]= 6.3618456410625559136428432181939059656629

In[5]:= N[(E^3 - 1)/3, 40]
Out[5]= 6.3618456410625559136428432181939059656629
```

## Implementation notes

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

**Attributes:** `NumericFunction`, `Protected`.

## References

- DLMF §13 — confluent hypergeometric functions (Kummer's 1F1).
- Source: [`src/special_functions/hypergeopfq.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/hypergeopfq.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_dsolve_m19_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m19_stress.c)
- Tests: [`tests/test_hypergeopfq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergeopfq.c)

## Notes & additional examples

### Notes

`Hypergeometric1F1[a, b, z]` is Kummer's confluent hypergeometric function, equal to `HypergeometricPFQ[{a}, {b}, z]`, and converges for all `z`. A non-positive integer `a` truncates the series to a polynomial (the Laguerre/Hermite family); otherwise the function evaluates numerically at machine, MPFR, and complex precision.
