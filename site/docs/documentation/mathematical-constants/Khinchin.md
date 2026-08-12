# Khinchin

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Khinchin`**

is Khinchin's constant K (also Khintchine's constant), with numerical value ~= 2.68545.

**`NumericQ[Khinchin] is True, and D[Khinchin, x] is 0. N[Khinchin, prec]`**

<details>
<summary>Notes</summary>

Khinchin's constant is the limiting geometric mean of the partial quotients in the continued-fraction expansion of almost every real number, given by the product over s \>= 1 of (1 + 1/(s (s + 2)))^Log2\[s\]. It is a mathematical constant: it has attributes Constant and Protected, evaluates it to any precision.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= N[Khinchin]
Out[1]= 2.68545
```

Khinchin's constant evaluates to arbitrary precision via its convergent product over partial quotients:

```mathematica
In[1]:= N[Khinchin, 60]
Out[1]= 2.685452001065306445309714835481795693820382293994462953051151
```

It is a true symbolic constant — `NumericQ` is `True` and its derivative vanishes:

```mathematica
In[1]:= NumericQ[Khinchin]
Out[1]= True

In[2]:= D[Khinchin, x]
Out[2]= 0
```

## Options & behaviour

The constant values for `GoldenAngle`, `Glaisher`, and `Khinchin` live in the
numeric constant table (`src/numeric.c`); their MPFR fillers compute
`GoldenAngle` from its closed form, and `Glaisher`/`Khinchin` from the series
above. Their `Constant`/`Protected` attributes are stamped in `numeric_init`.

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[Khinchin] = {Constant,
  Protected}`; the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[Khinchin]` is `True`
  and `D[Khinchin, x] = 0`.
- `N[Khinchin]` gives the machine value `2.68545`; `N[Khinchin, prec]` gives any
  precision, e.g.
  `N[Khinchin, 50] = 2.6854520010653064453097148354817956938203822939945`.

  Arbitrary precision uses the geometrically convergent zeta series
  `ln K · ln 2 = Σ_{n>=1} (ζ(2n) − 1)/n · Σ_{k=1}^{2n−1} (−1)^(k+1)/k`
  (the Bailey–Borwein–Crandall form). Verified to 250 digits.

**Attributes:** `Constant`, `Protected`.

## See also

[GoldenAngle](../../mathematical-constants/GoldenAngle/), [Glaisher](../../mathematical-constants/Glaisher/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)

## Notes & additional examples

### Notes

`Khinchin` is Khinchin's (Khintchine's) constant `K ~= 2.68545`, the limiting
geometric mean of the partial quotients in the continued-fraction expansion of
almost every real number: `K = Product[(1 + 1/(s (s + 2)))^Log2[s], {s, 1,
Infinity}]`. It carries the `Constant` and `Protected` attributes, so it stays
symbolic until `N[Khinchin, prec]` evaluates it to the requested precision.
