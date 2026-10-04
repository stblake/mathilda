# GoldenRatio

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`GoldenRatio`**

is the golden ratio phi = (1 + Sqrt\[5\])/2, with numerical value ~= 1.61803.

<details>
<summary>Notes</summary>

GoldenRatio is the positive root of x^2 == x + 1. It is a mathematical constant: it has attributes Constant and Protected, NumericQ\[GoldenRatio\] is True, and D\[GoldenRatio, x\] is 0. N\[GoldenRatio, prec\] evaluates it to any precision.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= N[GoldenRatio]
Out[1]= 1.61803

In[2]:= N[GoldenRatio, 40]
Out[2]= 1.6180339887498948482045868343656381177203

In[3]:= N[GoldenRatio^2 - GoldenRatio - 1, 40]
Out[3]= 1.9913648889155653462861004972096397569889e-59

In[4]:= FromContinuedFraction[{1, {1}}]
Out[4]= 1/2 (1 + Sqrt[5])

In[5]:= Round[N[(GoldenRatio^15 - (1 - GoldenRatio)^15)/Sqrt[5]]]
Out[5]= 610
```

## Implementation notes

**Definition.** `GoldenRatio` is the golden ratio φ = (1 + √5)/2 ≈ 1.6180339887, the
positive root of `x^2 == x + 1`. It is a constant *symbol*: row
`{ "GoldenRatio", 1.6180339887…, fill_mpfr_goldenratio }` of `kConstants[]` in
`src/numeric.c` supplies its value and the `numericalize_init` loop stamps it
`Constant | Protected`, so `Attributes[GoldenRatio]` is `{Constant, Protected}`.
Interned name `SYM_GoldenRatio` (`src/sym_names.c`), docstring in `src/info.c`.
`GoldenRatio` is a genuine symbol — it does *not* auto-expand to `(1+Sqrt[5])/2`
(`FullForm[GoldenRatio]` is `GoldenRatio`) — but it carries its defining relation
numerically, so `GoldenRatio^2 == GoldenRatio + 1` returns `True`. `D[GoldenRatio, x]`
is `0` from `Constant`.

**Representation & numeric value.** A bare `EXPR_SYMBOL`. `N[GoldenRatio]` returns the
tabled machine double; `N[GoldenRatio, k]` runs `fill_mpfr_goldenratio`, which forms
`(1 + sqrt(5))/2` in MPFR at `k + 20` guard bits and rounds down to the target
precision. Machine value 1.61803;
`N[GoldenRatio, 50] = 1.61803398874989484820458683436563811772030917980576`.

**Usage & limits.** `GoldenRatio` is the quadratic surd behind Fibonacci / Lucas
closed forms (Binet's formula) and the continued fraction `[1; 1, 1, …]`;
`NumericQ[GoldenRatio]` is `True` (whitelisted in `is_numeric_quantity`, `src/core.c`).
Because it is kept as a symbol rather than expanded to its radical form, algebraic
simplifications that need the surd must apply the defining relation (or `N`)
explicitly.

- Attributes `Constant`, `Protected`. `Attributes[GoldenRatio] = {Constant,
  Protected}`; the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[GoldenRatio]` is `True`
  and `D[GoldenRatio, x] = 0`.
- `N[GoldenRatio]` gives the machine value `1.61803`; `N[GoldenRatio, prec]`
  gives any precision, e.g.
  `N[GoldenRatio, 50] = 1.61803398874989484820458683436563811772030917980576`.

**Attributes:** `Constant`, `Protected`.

## References

- S. R. Finch, *Mathematical Constants* (Cambridge Univ. Press, 2003), §1.2.
- Source: [`src/numeric.c`](https://github.com/stblake/mathilda/blob/main/src/numeric.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)

## Notes & additional examples

### Notes

`GoldenRatio` is `phi = (1 + Sqrt[5])/2`, the positive root of `x^2 == x + 1`;
evaluating that polynomial at `phi` numerically returns `0.0`, confirming the
defining relation. It has the simplest possible continued fraction `[1; 1, 1, ...]`,
so `FromContinuedFraction[{1, {1}}]` recovers it exactly. Through Binet's formula
`Fibonacci[n] == (phi^n - (1 - phi)^n)/Sqrt[5]`, the closed form at `n = 15`
rounds to `610 = Fibonacci[15]`. It is a protected `Constant` (so
`D[GoldenRatio, x]` is `0`) evaluated to any precision by `N`.
