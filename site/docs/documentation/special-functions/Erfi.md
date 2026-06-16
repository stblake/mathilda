# Erfi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Erfi[z]
    gives the imaginary error function erfi(z) = -I Erf[I z]
    = (2/Sqrt[Pi]) Integral_0^z e^(t^2) dt.
Erfi[0] = 0, Erfi[Infinity] = Infinity, Erfi[I Infinity] = I. An entire
function, odd in z. Real and complex inputs evaluate numerically at
machine or arbitrary (MPFR) precision; D[Erfi[z], z] = (2/Sqrt[Pi]) E^(z^2).
Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Exact special values: `Erfi[0] = 0`, `Erfi[Infinity] = Infinity`,

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Erfi[0]
Out[1]= 0
```

```mathematica
In[1]:= Erfi[-x]
Out[1]= -Erfi[x]
```

```mathematica
In[1]:= N[Erfi[1], 30]
Out[1]= 1.650425758797542876025337729561
```

```mathematica
In[1]:= Series[Erfi[x], {x, 0, 7}]
Out[1]= 2/Sqrt[Pi] x + 2/3/Sqrt[Pi] x^3 + 1/5/Sqrt[Pi] x^5 + 1/21/Sqrt[Pi] x^7 + O[x]^8

In[2]:= D[Erfi[x], x]
Out[2]= (2 E^x^2)/Sqrt[Pi]
```

### Notes

`Erfi[z] = -I Erf[I z] = (2/Sqrt[Pi]) Integral_0^z e^(t^2) dt` is the imaginary
error function, an entire odd function with `Erfi[0] = 0`,
`Erfi[Infinity] = Infinity`, and `Erfi[I Infinity] = I`. Compared with `Erf`, the
sign of every Maclaurin coefficient is positive, reflecting the `+t^2` in the
integrand. Real and complex arguments evaluate numerically at machine or
arbitrary (MPFR) precision, the derivative is `(2/Sqrt[Pi]) E^(z^2)`, and `Erfi`
is `Listable`.
