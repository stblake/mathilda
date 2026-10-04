---
source: src/numeric.c
---
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
