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

In[2]:= N[Khinchin, 60]
Out[2]= 2.685452001065306445309714835481795693820382293994462953051151

In[3]:= NumericQ[Khinchin]
Out[3]= True

In[4]:= D[Khinchin, x]
Out[4]= 0
```

## Options & behaviour

The constant values for `GoldenAngle`, `Glaisher`, and `Khinchin` live in the
numeric constant table (`src/numeric.c`); their MPFR fillers compute
`GoldenAngle` from its closed form, and `Glaisher`/`Khinchin` from the series
above. Their `Constant`/`Protected` attributes are stamped in `numeric_init`.

## Implementation notes

**Definition.** `Khinchin` is Khinchin's constant *K* (also "Khintchine's constant")
≈ 2.6854520011: the geometric mean of the partial quotients in the continued-fraction
expansion of almost every real number, `Product_{s>=1} (1 + 1/(s (s + 2)))^Log2[s]`.
It is a constant *symbol*: row `{ "Khinchin", 2.6854520011…, fill_mpfr_khinchin }` of
`kConstants[]` in `src/numeric.c` supplies its value and the `numericalize_init` loop
stamps it `Constant | Protected`, so `Attributes[Khinchin]` is `{Constant, Protected}`.
Interned name `SYM_Khinchin` (`src/sym_names.c`), docstring in `src/info.c`. No
DownValues; `D[Khinchin, x] -> 0` follows from `Constant`.

**Representation & numeric value.** A bare `EXPR_SYMBOL`, kept exact under evaluation.
`N[Khinchin]` returns the tabled machine double. `N[Khinchin, k]` runs
`fill_mpfr_khinchin`, which uses the geometrically convergent **Bailey–Borwein–Crandall
zeta series** `Log[K] Log[2] = Sum_{n>=1} ((Zeta[2n] - 1)/n) Sum_{k=1}^{2n-1}
(-1)^(k+1)/k`: `Zeta[2n]-1 ~ 4^-n`, so ~`bits/2` terms (each via `mpfr_zeta_ui`, with a
running alternating-harmonic inner sum) reach the target; `K = Exp[S / Log[2]]` at
`bits + 64` guard precision. Machine value 2.68545;
`N[Khinchin, 50] = 2.68545200106530644530971483548179569382038229399446`.

**Usage & limits.** `Khinchin` is the metric constant of continued fractions;
`NumericQ[Khinchin]` is `True` (whitelisted in `is_numeric_quantity`, `src/core.c`).
Whether *K* is irrational is open. The implementation limit mirrors `Glaisher`'s: the
filler is a bespoke convergent series rather than a library constant, with a term
cutoff at `2^-(bits+32)`.

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

## References

**See also:** [GoldenAngle](../../mathematical-constants/GoldenAngle/), [Glaisher](../../mathematical-constants/Glaisher/)

- D. H. Bailey, J. M. Borwein and R. E. Crandall, *On the Khintchine constant*, Math. Comp. **66** (1997) 417-431.
- S. R. Finch, *Mathematical Constants* (Cambridge Univ. Press, 2003), §1.8.
- Source: [`src/numeric.c`](https://github.com/stblake/mathilda/blob/main/src/numeric.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)

## Notes & additional examples

### Notes

`Khinchin` is Khinchin's (Khintchine's) constant `K ~= 2.68545`, the limiting
geometric mean of the partial quotients in the continued-fraction expansion of
almost every real number: `K = Product[(1 + 1/(s (s + 2)))^Log2[s], {s, 1,
Infinity}]`. It carries the `Constant` and `Protected` attributes, so it stays
symbolic until `N[Khinchin, prec]` evaluates it to the requested precision.
