---
source: src/numeric.c
---
**Definition.** `Pi` is the circle constant π = 3.14159…, the ratio of a circle's
circumference to its diameter. It is a first-class constant *symbol*, not a function:
the `kConstants[]` table in `src/numeric.c` registers it with a machine value and an
MPFR filler, and the init loop at the bottom of `numericalize_init` stamps every row
of that table `Constant | Protected`, so `Attributes[Pi]` is `{Constant, Protected}`.
The interned name is `SYM_Pi` (`src/sym_names.c`); the docstring is in `src/info.c`.
`Pi` carries no DownValues of its own. The exact closed forms that *produce* it —
`Sin[Pi] -> 0`, `Cos[Pi] -> -1`, `ArcTan[1] -> Pi/4`, `Exp[I Pi] -> -1` — live in the
trigonometric and log-exp heads, and `D[Pi, x] -> 0` follows from the `Constant`
attribute.

**Representation & numeric value.** Internally `Pi` is a bare `EXPR_SYMBOL`; it
survives evaluation unevaluated and `FullForm[Pi]` is `Pi`. `N[Pi]` returns the
machine double `M_PI` (via `leaf_from_double`). `N[Pi, k]` goes through
`numericalize_symbol`, which calls the filler `fill_mpfr_pi` → MPFR's `mpfr_const_pi`
at the requested bit precision and returns an `EXPR_MPFR`. Machine value 3.14159;
`N[Pi, 50] = 3.1415926535897932384626433832795028841971693993751`.

**Usage & limits.** `Pi` drives the special-value reductions of the elementary and
inverse-trigonometric functions and appears throughout the Zeta / Gamma / series
results; `NumericQ[Pi]` is `True` (the symbol is whitelisted in `is_numeric_quantity`,
`src/core.c`). It is never reduced to a rational — being irrational, its only
"evaluation" is that its digits are materialised on demand by `N`.
