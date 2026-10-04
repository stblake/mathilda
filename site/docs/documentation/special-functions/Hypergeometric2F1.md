# Hypergeometric2F1

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Hypergeometric2F1[a, b, c, z]`**

is the Gauss hypergeometric 2F1, equal to HypergeometricPFQ\[{a, b}, {c}, z\].

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Hypergeometric2F1[1, 1, 2, z]
Out[1]= -Log[1 - z]/z
```

### Applications (4)

```mathematica
In[2]:= Hypergeometric2F1[1, 1, 2, z]
Out[2]= -Log[1 - z]/z

In[3]:= Hypergeometric2F1[-3, 1, 1, z]
Out[3]= 1 - 3 z + 3 z^2 - z^3

In[4]:= N[Hypergeometric2F1[1/2, 1/2, 3/2, 1/4]/2, 40]
Out[4]= 0.52359877559829887307710723054658381403285

In[5]:= N[ArcSin[1/2], 40]
Out[5]= 0.52359877559829887307710723054658381403285
```

## Implementation notes

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

**Attributes:** `NumericFunction`, `Protected`.

## References

- DLMF §15 — the Gauss hypergeometric function 2F1.
- Source: [`src/special_functions/hypergeopfq.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/hypergeopfq.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)

## Notes & additional examples

### Notes

`Hypergeometric2F1[a, b, c, z]` is the Gauss hypergeometric function `HypergeometricPFQ[{a, b}, {c}, z]`, convergent for `|z| < 1` (and by termination for non-positive integer `a` or `b`). Many elementary functions are special cases: `2F1[1, 1, 2, z] = -Log[1 - z]/z` and `z * 2F1[1/2, 1/2, 3/2, z^2] = ArcSin[z]`.
