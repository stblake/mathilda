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

**Attributes:** `NumericFunction`, `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_hypergeopfq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergeopfq.c)

## Notes & additional examples

### Notes

`Hypergeometric2F1[a, b, c, z]` is the Gauss hypergeometric function `HypergeometricPFQ[{a, b}, {c}, z]`, convergent for `|z| < 1` (and by termination for non-positive integer `a` or `b`). Many elementary functions are special cases: `2F1[1, 1, 2, z] = -Log[1 - z]/z` and `z * 2F1[1/2, 1/2, 3/2, z^2] = ArcSin[z]`.
