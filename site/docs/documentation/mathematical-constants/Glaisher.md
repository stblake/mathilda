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

In[2]:= N[Glaisher, 40]
Out[2]= 1.2824271291006226368753425688697917277676

In[3]:= NumericQ[Glaisher]
Out[3]= True

In[4]:= D[Glaisher, x]
Out[4]= 0
```

## Implementation notes

**Definition.** `Glaisher` is the Glaisher–Kinkelin constant *A* ≈ 1.2824271291, the
constant appearing in the asymptotics of the hyperfactorial and the Barnes *G*
function. It satisfies `Log[A] == 1/12 - Zeta'[-1]`, equivalently
`Log[A] = (gamma + Log[2 Pi])/12 - Zeta'[2]/(2 Pi^2)`. It is a constant *symbol*: row
`{ "Glaisher", 1.2824271291…, fill_mpfr_glaisher }` of `kConstants[]` in
`src/numeric.c` supplies its value and the `numericalize_init` loop stamps it
`Constant | Protected`, so `Attributes[Glaisher]` is `{Constant, Protected}`. Interned
name `SYM_Glaisher` (`src/sym_names.c`), docstring in `src/info.c`. No DownValues;
`D[Glaisher, x] -> 0` follows from `Constant`.

**Representation & numeric value.** A bare `EXPR_SYMBOL`, kept exact under evaluation.
`N[Glaisher]` returns the tabled machine double. `N[Glaisher, k]` runs the hand-rolled
`fill_mpfr_glaisher`: MPFR has no ζ′ primitive, so `Zeta'[2] = -Sum_{n>=1} Log[n]/n^2`
is evaluated by **Euler–Maclaurin summation** (head sum to `N ~ bits/6`, the integral
tail `(ln N + 1)/N`, the `f(N)/2` term, and a divergent asymptotic correction series
whose Bernoulli numbers come from `mpfr_zeta_ui`, truncated at its smallest term);
`Log[A]` is then assembled and exponentiated at `bits + 64` guard precision. Machine
value 1.28243; `N[Glaisher, 50] = 1.282427129100622636875342568869791727767688927325`.

**Usage & limits.** `Glaisher` is the closed form behind hyperfactorial / Barnes-G
asymptotics; `NumericQ[Glaisher]` is `True` (whitelisted in `is_numeric_quantity`,
`src/core.c`). The only implementation limit is that its filler is a bespoke
Euler–Maclaurin routine (not a library constant): the guard bits and smallest-term
truncation are sized so the asymptotic floor lies below any requested precision.

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

## References

**See also:** [Zeta](../../special-functions/Zeta/)

- S. R. Finch, *Mathematical Constants* (Cambridge Univ. Press, 2003), §2.15.
- Source: [`src/numeric.c`](https://github.com/stblake/mathilda/blob/main/src/numeric.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)

## Notes & additional examples

### Notes

`Glaisher` is the Glaisher-Kinkelin constant `A`, defined by
`Log[A] == 1/12 - Zeta'[-1]` and appearing in the asymptotics of the
hyperfactorial and in many `Zeta`-derivative identities. It is held symbolic
(attributes `Constant` and `Protected`, so `D[Glaisher, x]` is `0` and
`NumericQ` is `True`) until `N` forces a value; `N[Glaisher, 40]` returns it
to 40 digits via its MPFR series.
