# EulerGamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EulerGamma`**

is Euler's constant gamma, with numerical value ~= 0.5772156649.

**`HarmonicNumber[n] - Log[n] as n -> Infinity. It is a mathematical`**

<details>
<summary>Notes</summary>

EulerGamma is the Euler-Mascheroni constant, the limit of constant: it has attributes Constant and Protected, NumericQ\[EulerGamma\] is True, and D\[EulerGamma, x\] is 0. N\[EulerGamma, prec\] evaluates it to any precision.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= N[EulerGamma]
Out[1]= 0.577216

In[2]:= N[EulerGamma, 50]
Out[2]= 0.577215664901532860606512090082402431042159335939923

In[3]:= D[EulerGamma, x]
Out[3]= 0

In[4]:= NumericQ[EulerGamma]
Out[4]= True

In[5]:= N[Gamma[1/2] + EulerGamma, 30]
Out[5]= 2.349669515807048887904679573425
```

## Options & behaviour

The constant value lives in the central numeric constant table (`src/numeric.c`)
the symbol's identity (attributes) is registered in `src/special_functions/eulergamma.c`.

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[EulerGamma] = {Constant,
  Protected}`; the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[EulerGamma]` is
  `True` and `D[EulerGamma, x] = 0`.
- `N[EulerGamma]` gives the machine value `0.577216`; `N[EulerGamma, prec]`
  gives any precision (MPFR `mpfr_const_euler`), e.g.
  `N[EulerGamma, 50] = 0.57721566490153286060651209008240243104215933593992`.
- Participates in exact numeric work, e.g.
  `Round[1/EulerGamma^100] = 734833795660954410469466`, and digit/continued-
  fraction extraction, e.g. `RealDigits[EulerGamma, 10, 50, -10^4]` returns
  decimal digits 10000–10049.

**Attributes:** `Constant`, `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)

## Notes & additional examples

### Notes

`EulerGamma` is the Euler-Mascheroni constant, the limit of
`HarmonicNumber[n] - Log[n]`. It carries the `Constant` and `Protected`
attributes, so it is `NumericQ` and differentiates to `0`, yet stays exact
until `N[EulerGamma, prec]` evaluates it to any requested precision via MPFR.
