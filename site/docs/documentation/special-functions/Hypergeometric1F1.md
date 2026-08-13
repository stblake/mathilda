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
Out[4]= 6.3618456410625559136428432181939059656621

In[5]:= N[(E^3 - 1)/3, 40]
Out[5]= 6.3618456410625559136428432181939059656621
```

## Implementation notes

**Attributes:** `NumericFunction`, `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_hypergeopfq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergeopfq.c)
- Tests: [`tests/test_integrate_ramanujan.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_ramanujan.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)

## Notes & additional examples

### Notes

`Hypergeometric1F1[a, b, z]` is Kummer's confluent hypergeometric function, equal to `HypergeometricPFQ[{a}, {b}, z]`, and converges for all `z`. A non-positive integer `a` truncates the series to a polynomial (the Laguerre/Hermite family); otherwise the function evaluates numerically at machine, MPFR, and complex precision.
