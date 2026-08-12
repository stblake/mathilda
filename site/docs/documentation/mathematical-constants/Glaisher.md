# Glaisher

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Glaisher`**

is the Glaisher-Kinkelin constant A, with numerical value ~= 1.28243.

**`D[Glaisher, x] is 0. N[Glaisher, prec] evaluates it to any precision.`**

<details>
<summary>Notes</summary>

Glaisher's constant satisfies Log\[A\] == 1/12 - Zeta'\[-1\], where Zeta is the Riemann zeta function. It is a mathematical constant: it has attributes Constant and Protected, NumericQ\[Glaisher\] is True, and

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= N[Glaisher]
Out[1]= 1.28243
```

```mathematica
In[1]:= N[Glaisher, 40]
Out[1]= 1.2824271291006226368753425688697917277676
```

```mathematica
In[1]:= NumericQ[Glaisher]
Out[1]= True
```

```mathematica
In[1]:= D[Glaisher, x]
Out[1]= 0
```

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[Glaisher] = {Constant,
  Protected}`; the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[Glaisher]` is `True`
  and `D[Glaisher, x] = 0`.
- `N[Glaisher]` gives the machine value `1.28243`; `N[Glaisher, prec]` gives any
  precision, e.g.
  `N[Glaisher, 50] = 1.282427129100622636875342568869791727767688927325`.

  Arbitrary precision is computed from `ln A = (γ + ln(2π))/12 − ζ'(2)/(2π²)`,
  with `ζ'(2)` evaluated by Euler–Maclaurin summation of `−Σ ln(n)/n²` (the
  even Bernoulli numbers obtained from `Zeta[2k]`). Verified to 250 digits.

**Attributes:** `Constant`, `Protected`.

## See also

[Zeta](../../special-functions/Zeta/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)

## Notes & additional examples

### Notes

`Glaisher` is the Glaisher-Kinkelin constant `A`, defined by
`Log[A] == 1/12 - Zeta'[-1]` and appearing in the asymptotics of the
hyperfactorial and in many `Zeta`-derivative identities. It is held symbolic
(attributes `Constant` and `Protected`, so `D[Glaisher, x]` is `0` and
`NumericQ` is `True`) until `N` forces a value; `N[Glaisher, 40]` returns it
to 40 digits via its MPFR series.
