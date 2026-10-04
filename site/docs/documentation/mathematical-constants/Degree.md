# Degree

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Degree`**

gives the number of radians in one degree, with numerical value Pi/180 (~= 0.0174533).

<details>
<summary>Notes</summary>

Multiply by Degree to convert degrees to radians, so 30 Degree is 30 degrees. It is a mathematical constant: it has attributes Constant and Protected, NumericQ\[Degree\] is True, and D\[Degree, x\] is 0. N\[Degree, prec\] evaluates it to any precision.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= N[Degree]
Out[1]= 0.0174533

In[2]:= Sin[30 Degree] // N
Out[2]= 0.5

In[3]:= N[Tan[60 Degree], 40]
Out[3]= 1.7320508075688772935274463415058723669427

In[4]:= N[Degree, 40]
Out[4]= 0.017453292519943295769236907684886127134428
```

## Options & behaviour

The constant values for `Pi`, `E`, `Catalan`, `GoldenRatio`, and `Degree` all
live in the central numeric constant table (`src/numeric.c`); their identities
(attributes `Constant`, `Protected`) are stamped in `numeric_init`.

## Implementation notes

**Definition.** `Degree` is the number of radians in one degree, i.e. the constant
π/180 ≈ 0.0174533. Multiplying an angle in degrees by `Degree` converts it to radians,
so `30 Degree` *represents* 30°. It is a constant *symbol*: row
`{ "Degree", M_PI/180.0, fill_mpfr_degree }` of the `kConstants[]` table in
`src/numeric.c` supplies its value and the `numericalize_init` loop stamps it
`Constant | Protected`, so `Attributes[Degree]` is `{Constant, Protected}`. Interned
name `SYM_Degree` (`src/sym_names.c`), docstring in `src/info.c`. `D[Degree, x] -> 0`
follows from `Constant`. `Degree` participates in the trigonometric closed forms
through its numeric relation to `Pi` rather than through a rewrite rule of its own:
`180 Degree == Pi` returns `True` and `Sin[90 Degree] -> 1`, while the product
`30 Degree` itself is left as the unevaluated `Times[30, Degree]`.

**Representation & numeric value.** `Degree` is a bare `EXPR_SYMBOL`. `N[Degree]`
returns the machine double `M_PI/180`; `N[Degree, k]` runs `fill_mpfr_degree`, which
computes `mpfr_const_pi` at `k + 10` guard bits, divides by `180`, and rounds down to
the target precision. Machine value 0.0174533;
`N[Degree, 50] = 0.0174532925199432957692369076848861271344287188854173`.

**Usage & limits.** `Degree` is the degrees-to-radians scale factor; `NumericQ[Degree]`
is `True` (whitelisted in `is_numeric_quantity`, `src/core.c`). A symbolic
`Cos[30 Degree]` is left unevaluated, but `N[Cos[30 Degree]]` evaluates to `0.866025`.
The only limit is the one shared by the `Pi`-derived constants: it is irrational and
its digits are produced on demand by `N`.

- Attributes `Constant`, `Protected`. `Attributes[Degree] = {Constant,
  Protected}`; the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[Degree]` is `True` and
  `D[Degree, x] = 0`.
- Used in arguments of trigonometric functions to express angles in degrees,
  e.g. `30 Degree` is `π/6`; the trig value evaluates numerically under `N`,
  e.g. `N[Cos[30 Degree]] = 0.866025` (the exact symbolic form
  `Cos[30 Degree]` is left unevaluated).
- `N[Degree]` gives the machine value `0.0174533`; `N[Degree, prec]` gives any
  precision, e.g.
  `N[Degree, 50] = 0.0174532925199432957692369076848861271344287188854173`.

**Attributes:** `Constant`, `Protected`.

## References

**See also:** [N](../../arithmetic/N/), [Pi](../../mathematical-constants/Pi/), [E](../../mathematical-constants/E/), [Catalan](../../mathematical-constants/Catalan/), [GoldenRatio](../../mathematical-constants/GoldenRatio/)

- Source: [`src/numeric.c`](https://github.com/stblake/mathilda/blob/main/src/numeric.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)

## Notes & additional examples

### Notes

`Degree` is the mathematical constant `Pi/180`, the number of radians in one degree. Multiply an angle by it to convert degrees to radians, so `30 Degree` is the radian measure of 30 degrees and `Sin[30 Degree]` numericalises to `0.5`. It carries the `Constant` and `Protected` attributes (`D[Degree, x]` is `0`), is recognised by `NumericQ`, and evaluates to arbitrary precision under `N[Degree, prec]` — the last example agrees with `N[Pi/180, 40]` digit for digit.
