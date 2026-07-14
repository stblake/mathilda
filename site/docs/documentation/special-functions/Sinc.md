# Sinc

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Sinc[z]
    gives the cardinal sine Sin[z]/z, with Sinc[0] = 1.
An entire, even function. Sinc[+-Infinity] = 0. Real and complex inputs
evaluate numerically at machine or arbitrary (MPFR) precision;
D[Sinc[z], z] = Cos[z]/z - Sin[z]/z^2. Listable.
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
In[1]:= Sinc[0]
Out[1]= 1
```

```mathematica
In[1]:= Sinc[2.]
Out[1]= 0.454649
```

```mathematica
In[1]:= N[Sinc[2], 45]
Out[1]= 0.454648713412840847698009932955872421351127485
```

```mathematica
In[1]:= Sinc[1. + I]
Out[1]= 0.966711 - 0.331747 I
```

```mathematica
In[1]:= D[Sinc[x], x]
Out[1]= Cos[x]/x - Sin[x]/x^2
```

```mathematica
In[1]:= Series[Sinc[x], {x, 0, 6}]
Out[1]= 1 - x^2/6 + x^4/120 - x^6/5040 + O[x]^7
```

### Notes

`Sinc[z]` is the cardinal sine `Sin[z]/z`, with the removable singularity at the
origin filled in as `Sinc[0] = 1`. It is entire and even, and `Sinc[±Infinity] = 0`.
It appears as the derivative of the sine integral: `D[SinIntegral[z], z] = Sinc[z]`.
Numeric evaluation is at machine or arbitrary (MPFR) precision for both real and
complex arguments. Listable. See also [`SinIntegral`](SinIntegral.md).
