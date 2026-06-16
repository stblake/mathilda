# StieltjesGamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StieltjesGamma[n]
    gives the n-th Stieltjes constant gamma_n, the Laurent coefficients of
Zeta about s = 1. StieltjesGamma[0] is EulerGamma; higher constants are
inert (they stay symbolic) and appear in Series expansions of Zeta. Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StieltjesGamma[0]
Out[1]= EulerGamma
```

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)

## Notes & additional examples

### Worked examples

The zeroth Stieltjes constant is the Euler–Mascheroni constant:

```mathematica
In[1]:= StieltjesGamma[0]
Out[1]= EulerGamma
```

Higher constants are inert and stay symbolic:

```mathematica
In[1]:= StieltjesGamma[3]
Out[1]= StieltjesGamma[3]
```

Numericalizing the zeroth constant recovers γ to 40 digits:

```mathematica
In[1]:= N[StieltjesGamma[0], 40]
Out[1]= 0.57721566490153286060651209008240243104214
```

The constants are exactly the Laurent coefficients of `Zeta` about `s = 1` —
expanding the series exhibits them in their defining role:

```mathematica
In[1]:= Series[Zeta[s], {s, 1, 2}]
Out[1]= 1/(s - 1) + EulerGamma + -StieltjesGamma[1] (s - 1) + 1/2 StieltjesGamma[2] (s - 1)^2 + O[s - 1]^3
```

### Notes

`StieltjesGamma[n]` denotes the `n`-th Stieltjes constant γ_n, defined by the
Laurent expansion `Zeta[s] = 1/(s - 1) + Sum[(-1)^n/n! γ_n (s - 1)^n]` about the
pole at `s = 1`. `StieltjesGamma[0]` is `EulerGamma`; the higher constants are
inert symbols that appear, with the correct `(-1)^n/n!` factors, in the `Series`
expansion of `Zeta` at `s = 1`. It is `Listable`.
