# Erf

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Erf[z]
    gives the error function erf(z) = (2/Sqrt[Pi]) Integral_0^z e^(-t^2) dt.
Erf[z0, z1]
    gives the generalized error function erf(z1) - erf(z0).
Erf[0] = 0, Erf[Infinity] = 1, Erf[-Infinity] = -1. An entire function,
odd in z. Real and complex inputs evaluate numerically at machine or
arbitrary (MPFR) precision; D[Erf[z], z] = (2/Sqrt[Pi]) E^(-z^2). Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Exact special values: `Erf[0] = 0`, `Erf[Infinity] = 1`,

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Erf[0]
Out[1]= 0
```

```mathematica
In[1]:= Erf[Infinity]
Out[1]= 1

In[2]:= Erf[-z]
Out[2]= -Erf[z]
```

```mathematica
In[1]:= N[Erf[1], 40]
Out[1]= 0.84270079294971486934122063508260925929605
```

```mathematica
In[1]:= N[Erf[1 + I], 20]
Out[1]= 1.31615128169794764488 + 0.190453469237834686284*I
```

```mathematica
In[1]:= Series[Erf[x], {x, 0, 7}]
Out[1]= 2/Sqrt[Pi] x + -2/3/Sqrt[Pi] x^3 + 1/5/Sqrt[Pi] x^5 + -1/21/Sqrt[Pi] x^7 + O[x]^8

In[2]:= D[Erf[x^2], x]
Out[2]= (4 x E^(-x^4))/Sqrt[Pi]
```

### Notes

`Erf[z]` is the error function `(2/Sqrt[Pi]) Integral_0^z e^(-t^2) dt`, an entire
odd function with the exact values `Erf[0] = 0` and `Erf[±Infinity] = ±1`. Real
and complex arguments evaluate numerically at machine or arbitrary (MPFR)
precision — the complex path uses a DLMF series so `N[Erf[1 + I], 20]` is correct
to the requested digits. The Maclaurin series and the chain-rule derivative
`D[Erf[z], z] = (2/Sqrt[Pi]) E^(-z^2)` are both built in, and `Erf` is `Listable`.
