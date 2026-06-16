# ExpIntegralEi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ExpIntegralEi[z]
    gives the exponential integral Ei(z), the principal value of
    -Integral_{-z}^Infinity e^-t/t dt, with a branch cut on (-Infinity, 0).
ExpIntegralEi[0] = -Infinity, ExpIntegralEi[Infinity] = Infinity,
ExpIntegralEi[-Infinity] = 0, ExpIntegralEi[+-I Infinity] = +-I Pi. Real and
complex inputs evaluate numerically at machine or arbitrary (MPFR) precision;
D[ExpIntegralEi[z], z] = E^z/z. Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Exact special values: `ExpIntegralEi[0] = -Infinity`,

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= ExpIntegralEi[0]
Out[1]= -Infinity
```

```mathematica
In[1]:= D[ExpIntegralEi[z], z]
Out[1]= E^z/z
```

```mathematica
In[1]:= N[ExpIntegralEi[1], 40]
Out[1]= 1.8951178163559367554665209343316342690171
```

```mathematica
In[1]:= N[ExpIntegralEi[I], 30]
Out[1]= 0.3374039229009681346626462038893 + 2.516879397162079634172675005462*I
```

### Notes

`ExpIntegralEi[z]` is the exponential integral `Ei(z)`, with a branch cut on
`(-Infinity, 0)` and derivative `E^z/z`. On the imaginary axis it ties to the
cosine/sine integrals via `Ei(I) = Ci(1) + I (Pi/2 + Si(1))`. Real and complex
arguments evaluate at machine or arbitrary (MPFR) precision. Listable.
