# CosIntegral

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
CosIntegral[z]
    gives the cosine integral Ci(z) = -Integral_z^Infinity Cos[t]/t dt.
Has a logarithmic singularity at 0 and a branch cut on (-Infinity, 0].
CosIntegral[0] = -Infinity, CosIntegral[Infinity] = 0,
CosIntegral[-Infinity] = I Pi, CosIntegral[+-I Infinity] = Infinity.
Real and complex inputs evaluate numerically at machine or arbitrary (MPFR)
precision; D[CosIntegral[z], z] = Cos[z]/z. Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= CosIntegral[2.8]
Out[1]= 0.186488
```

```mathematica
In[1]:= N[CosIntegral[2], 50]
Out[1]= 0.42298082877486499569856515319825589413573775630619
```

```mathematica
In[1]:= CosIntegral[{-Infinity, Infinity, -I Infinity, I Infinity}]
Out[1]= {I Pi, 0, Infinity, Infinity}
```

```mathematica
In[1]:= CosIntegral[-2.]
Out[1]= 0.422981 + 3.14159 I
```

```mathematica
In[1]:= CosIntegral[3. I]
Out[1]= 4.96039 + 1.5708 I
```

```mathematica
In[1]:= D[CosIntegral[x], x]
Out[1]= Cos[x]/x
```

```mathematica
In[1]:= Series[CosIntegral[x], {x, 0, 6}]
Out[1]= EulerGamma + Log[x] - 1/4 x^2 + 1/96 x^4 - 1/4320 x^6 + O[x]^7
```

```mathematica
In[1]:= Normal[Series[CosIntegral[x], {x, Infinity, 3}]]
Out[1]= -Cos[x]/x^2 + Sin[x] (1/x - 2/x^3)
```

### Notes

`CosIntegral[z]` is the cosine integral `Ci(z) = -Integral_z^Infinity Cos[t]/t dt`.
Unlike its sibling [`SinIntegral`](SinIntegral.md) — which is entire and odd — `Ci`
has a logarithmic singularity at the origin (`CosIntegral[0] = -Infinity`) and a
branch cut running along the negative real axis `(-Infinity, 0]`. On the cut it
takes the from-above value, so for a negative real `x` the result is complex:
`Ci(x) = Ci(|x|) + I Pi`, matching `CosIntegral[-Infinity] = I Pi`. On the imaginary
axis `Ci(I y) = Chi(y) + I Pi/2` in terms of the hyperbolic cosine integral. Its
derivative is `Cos[z]/z`. Numeric evaluation uses a convergent series near the
origin and a trig-prefactored asymptotic expansion for large `|z|`, at machine or
arbitrary (MPFR) precision. Listable.
