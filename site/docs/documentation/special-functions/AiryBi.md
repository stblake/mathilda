# AiryBi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
AiryBi[z]
    gives the Airy function Bi(z), the solution of y'' = z y that grows
    exponentially as z -> +Infinity.
AiryBi[0] = 1/(3^(1/6) Gamma[2/3]), AiryBi[Infinity] = Infinity,
AiryBi[-Infinity] = 0. An entire function of z. Real and complex inputs
evaluate numerically at machine or arbitrary (MPFR) precision;
D[AiryBi[z], z] = AiryBiPrime[z]. Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= AiryBi[0]
Out[1]= 1/(3^(1/6) Gamma[2/3])

In[2]:= AiryBi[1.8]
Out[2]= 2.59587
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= AiryBi[0]
Out[1]= 1/(3^(1/6) Gamma[2/3])
```

```mathematica
In[1]:= N[AiryBi[0], 40]
Out[1]= 0.6149266274460007351509223690936135535947
```

```mathematica
In[1]:= D[AiryBi[z], z]
Out[1]= AiryBiPrime[z]
```

```mathematica
In[1]:= N[AiryBi[2.0 + 1.0 I], 20]
Out[1]= 0.778230383757041677129 + 2.50509630006410244363*I
```

### Notes

`AiryBi[z]` is the dominant solution of the Airy equation `y'' == z y`, growing
exponentially as `z -> +Infinity` while `AiryBi[-Infinity] == 0`. Its exact
value at the origin is `1/(3^(1/6) Gamma[2/3])`, and `D[AiryBi[z], z]` returns
`AiryBiPrime[z]`. Complex arguments are evaluated to the requested MPFR
precision; `AiryBi` is Listable.
