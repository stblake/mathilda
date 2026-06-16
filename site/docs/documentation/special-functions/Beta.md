# Beta

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Beta[a, b]
    is the Euler beta function B(a, b) = Gamma(a) Gamma(b) / Gamma(a+b).
Beta[z, a, b]
    is the incomplete beta function Integral_0^z t^(a-1) (1-t)^(b-1) dt.
Beta[z0, z1, a, b]
    is the generalized incomplete beta Beta[z1, a, b] - Beta[z0, a, b].
Exact for rational arguments (a positive integer gives a rational via
Pochhammer); non-positive integer poles give ComplexInfinity. Machine
and arbitrary-precision (MPFR) real and complex inputs evaluate
numerically. The incomplete form reduces through Hypergeometric2F1.
Listable.
```

## Examples

_No verified examples yet for this function._

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
In[1]:= Beta[3, 5]
Out[1]= 1/105
```

The central value is exactly `Pi`, and rational orders fold into Gamma quotients:

```mathematica
In[1]:= Beta[1/2, 1/2]
Out[1]= Pi

In[2]:= Beta[1/3, 1/3]
Out[2]= Gamma[1/3]^2/Gamma[2/3]
```

Positive-integer orders give the reciprocal binomial relation `1/B(7, 3) = 9 C(8, 2)`:

```mathematica
In[1]:= Beta[7, 3]
Out[1]= 1/252

In[2]:= 1/Beta[7, 3] - 9 Binomial[8, 2]
Out[2]= 0
```

The three-argument incomplete beta, and arbitrary-precision numerics:

```mathematica
In[1]:= Beta[5, 2, 3]
Out[1]= 1025/12

In[2]:= N[Beta[2.5, 3.5], 30]
Out[2]= 0.036815538909255388078101134397
```

### Notes

`Beta[a, b] = Gamma[a] Gamma[b]/Gamma[a+b]` is the Euler beta function. `Beta[z, a, b]` is the incomplete beta integral, and `Beta[z0, z1, a, b]` the generalized incomplete form. Exact for rational arguments via Pochhammer; non-positive integer poles give `ComplexInfinity`; the incomplete form reduces through `Hypergeometric2F1`. Real and complex inputs evaluate at machine or MPFR precision. Listable.
