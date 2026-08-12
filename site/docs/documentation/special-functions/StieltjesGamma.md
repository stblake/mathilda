# StieltjesGamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StieltjesGamma[n]`**

gives the n-th Stieltjes constant gamma\_n, the Laurent coefficients of

<details>
<summary>Notes</summary>

Zeta about s = 1. StieltjesGamma\[0\] is EulerGamma; higher constants are inert (they stay symbolic) and appear in Series expansions of Zeta. Listable.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= StieltjesGamma[0]
Out[1]= EulerGamma
```

### Applications (4)

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

## Algorithm

Mathilda -- StieltjesGamma, the Stieltjes constants gamma_n.

```text
  StieltjesGamma[n] = gamma_n, the coefficients of the Laurent expansion
  of the Riemann zeta function about s = 1:

    zeta(s) = 1/(s-1) + Sum_{n>=0} ((-1)^n / n!) gamma_n (s-1)^n.
```

gamma_0 is the Euler-Mascheroni constant EulerGamma. The higher constants have no elementary closed form, so StieltjesGamma is inert: it stays symbolic, except for the single reduction StieltjesGamma[0] -> EulerGamma. It is the natural output of Series[Zeta[x], {x, 1, n}] (and the Taylor expansion of Zeta about 0).

Like LogGamma (see src/polygamma.c), this module owns only the symbol's identity and the n = 0 reduction; all generic symbolic behaviour comes from the evaluator.

Attributes: Listable, Protected.

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## See also

[Zeta](../../special-functions/Zeta/), [EulerGamma](../../mathematical-constants/EulerGamma/), [Series](../../power-series/Series/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/tests/test_flint_bridge.c)
- Tests: [`tests/test_residue.c`](https://github.com/stblake/mathilda/blob/main/tests/test_residue.c)
- Tests: [`tests/test_stieltjesgamma.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stieltjesgamma.c)
- Tests: [`tests/test_zeta.c`](https://github.com/stblake/mathilda/blob/main/tests/test_zeta.c)

## Notes & additional examples

### Notes

`StieltjesGamma[n]` denotes the `n`-th Stieltjes constant γ_n, defined by the
Laurent expansion `Zeta[s] = 1/(s - 1) + Sum[(-1)^n/n! γ_n (s - 1)^n]` about the
pole at `s = 1`. `StieltjesGamma[0]` is `EulerGamma`; the higher constants are
inert symbols that appear, with the correct `(-1)^n/n!` factors, in the `Series`
expansion of `Zeta` at `s = 1`. It is `Listable`.
