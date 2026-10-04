---
source: src/numeric.c
---
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
