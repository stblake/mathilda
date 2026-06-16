# AiryAi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
AiryAi[z]
    gives the Airy function Ai(z), the solution of y'' = z y that decays
    as z -> +Infinity.
AiryAi[0] = 1/(3^(2/3) Gamma[2/3]), AiryAi[+-Infinity] = 0. An entire
function of z. Real and complex inputs evaluate numerically at machine
or arbitrary (MPFR) precision; D[AiryAi[z], z] = AiryAiPrime[z]. Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= AiryAi[0]
Out[1]= 1/(3^(2/3) Gamma[2/3])

In[2]:= AiryAi[1.8]
Out[2]= 0.0470362
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
In[1]:= AiryAi[0]
Out[1]= 1/(3^(2/3) Gamma[2/3])
```

```mathematica
In[1]:= N[AiryAi[1], 40]
Out[1]= 0.13529241631288141552414742351546630617494
```

```mathematica
In[1]:= D[AiryAi[z], {z, 2}]
Out[1]= z AiryAi[z]
```

```mathematica
In[1]:= Series[AiryAi[z], {z, 0, 4}]
Out[1]= 1/(3^(2/3) Gamma[2/3]) + -1/(3^(1/3) Gamma[1/3]) z + 1/6/(3^(2/3) Gamma[2/3]) z^3 + -1/12/(3^(1/3) Gamma[1/3]) z^4 + O[z]^5
```

```mathematica
In[1]:= N[AiryAi[2] AiryBiPrime[2] - AiryAiPrime[2] AiryBi[2], 30]
Out[1]= 0.3183098861837906715377675267449
```

### Notes

`AiryAi[z]` is the recessive solution of the Airy equation `y'' == z y` decaying
as `z -> +Infinity`; the second-derivative identity `D[AiryAi[z], {z, 2}] == z AiryAi[z]`
makes this explicit. The exact origin value is `1/(3^(2/3) Gamma[2/3])`, and the
Maclaurin series shows the characteristic missing `z^2` term (every third
coefficient vanishes). The last example is the Airy Wronskian
`Ai(z) Bi'(z) - Ai'(z) Bi(z) == 1/Pi`, recovered numerically as
`0.318309886...`. Real and complex arguments evaluate at machine or MPFR
precision; `AiryAi` is Listable.
