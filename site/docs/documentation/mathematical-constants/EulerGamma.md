# EulerGamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
EulerGamma
    is Euler's constant gamma, with numerical value ~= 0.5772156649.
EulerGamma is the Euler-Mascheroni constant, the limit of
HarmonicNumber[n] - Log[n] as n -> Infinity. It is a mathematical
constant: it has attributes Constant and Protected, NumericQ[EulerGamma]
is True, and D[EulerGamma, x] is 0. N[EulerGamma, prec] evaluates it to
any precision.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[EulerGamma] = {Constant,

**Attributes:** `Constant`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= N[EulerGamma]
Out[1]= 0.577216
```

```mathematica
In[1]:= N[EulerGamma, 50]
Out[1]= 0.577215664901532860606512090082402431042159335939923
```

```mathematica
In[1]:= D[EulerGamma, x]
Out[1]= 0

In[2]:= NumericQ[EulerGamma]
Out[2]= True
```

```mathematica
In[1]:= N[Gamma[1/2] + EulerGamma, 30]
Out[1]= 2.349669515807048887904679573425
```

### Notes

`EulerGamma` is the Euler-Mascheroni constant, the limit of
`HarmonicNumber[n] - Log[n]`. It carries the `Constant` and `Protected`
attributes, so it is `NumericQ` and differentiates to `0`, yet stays exact
until `N[EulerGamma, prec]` evaluates it to any requested precision via MPFR.
