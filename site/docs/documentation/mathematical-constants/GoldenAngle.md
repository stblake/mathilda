# GoldenAngle

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GoldenAngle`**

is the golden angle (3 - Sqrt\[5\]) Pi = 2 Pi / GoldenRatio^2, with numerical value ~= 2.39996 radians (~= 137.5 degrees).

**`N[GoldenAngle, prec] evaluates it to any precision.`**

<details>
<summary>Notes</summary>

GoldenAngle is a mathematical constant: it has attributes Constant and Protected, NumericQ\[GoldenAngle\] is True, and D\[GoldenAngle, x\] is 0.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= N[GoldenAngle]
Out[1]= 2.39996

In[2]:= N[GoldenAngle, 40]
Out[2]= 2.399963229728653322231555506633613853125

In[3]:= N[GoldenAngle/Degree, 30]
Out[3]= 137.5077640500378546463487396284

In[4]:= N[GoldenAngle - (3 - Sqrt[5]) Pi, 40]
Out[4]= 3.9827297778311306925722009944192795139778e-59
```

## Implementation notes

**Definition.** `GoldenAngle` is the golden angle, the closed form `(3 - Sqrt[5]) Pi`
= `2 Pi / GoldenRatio^2` ≈ 2.39996 radians (≈ 137.5°) — the smaller of the two arcs
that divide a circle in the golden ratio, the spacing that appears in phyllotaxis. It
is a constant *symbol*: row `{ "GoldenAngle", 2.3999632297…, fill_mpfr_goldenangle }`
of `kConstants[]` in `src/numeric.c` supplies its value and the `numericalize_init`
loop stamps it `Constant | Protected`, so `Attributes[GoldenAngle]` is
`{Constant, Protected}`. Interned name `SYM_GoldenAngle` (`src/sym_names.c`), docstring
in `src/info.c`. It carries its defining relation numerically —
`GoldenAngle == 2 Pi / GoldenRatio^2` returns `True` — and `D[GoldenAngle, x] -> 0`
follows from `Constant`.

**Representation & numeric value.** A bare `EXPR_SYMBOL` that stays unevaluated
(`GoldenAngle` prints as itself). `N[GoldenAngle]` returns the tabled machine double;
`N[GoldenAngle, k]` runs `fill_mpfr_goldenangle`, which forms `(3 - sqrt(5)) * pi` in
MPFR at `k + 20` guard bits (using `mpfr_const_pi`) and rounds down to the target
precision. Machine value 2.39996;
`N[GoldenAngle, 50] = 2.39996322972865332223155550663361385312499901105812`.

**Usage & limits.** `GoldenAngle` is primarily a geometric constant; `NumericQ[GoldenAngle]`
is `True` (whitelisted in `is_numeric_quantity`, `src/core.c`). Like the other
`Pi`-derived constants it is irrational and is kept symbolic until `N` materialises its
digits.

- Attributes `Constant`, `Protected`. `Attributes[GoldenAngle] = {Constant,
  Protected}`; the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[GoldenAngle]` is `True`
  and `D[GoldenAngle, x] = 0`.
- `N[GoldenAngle]` gives the machine value `2.39996`; `N[GoldenAngle, prec]`
  gives any precision (computed in MPFR from the closed form `(3 - Sqrt[5])
  Pi`), e.g.
  `N[GoldenAngle, 50] = 2.3999632297286533222315555066336138531249990110581`.

**Attributes:** `Constant`, `Protected`.

## References

- Source: [`src/numeric.c`](https://github.com/stblake/mathilda/blob/main/src/numeric.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)

## Notes & additional examples

### Notes

`GoldenAngle` is `(3 - Sqrt[5]) Pi = 2 Pi / GoldenRatio^2`, the angle that
divides a full turn in the golden ratio — about `137.5` degrees, the divergence
angle that governs optimal phyllotactic spiral packing in plants. Dividing by
`Degree` shows the familiar `137.5...`, and subtracting the closed form
`(3 - Sqrt[5]) Pi` numerically returns `0.0`. It is a protected `Constant`
(so `D[GoldenAngle, x]` is `0`) that `N` evaluates to any precision.
