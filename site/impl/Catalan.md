---
source: src/numeric.c
references:
  - "S. R. Finch, *Mathematical Constants* (Cambridge Univ. Press, 2003), §1.7."
---
**Definition.** `Catalan` is Catalan's constant *G* ≈ 0.9159655942, the alternating
sum `Sum_{k>=0} (-1)^k (2 k + 1)^-2` (equivalently `Beta[2]`, the Dirichlet beta
function at 2). It is a constant *symbol*: row
`{ "Catalan", 0.9159655941…, fill_mpfr_catalan }` of `kConstants[]` in
`src/numeric.c` supplies its value and the `numericalize_init` loop stamps it
`Constant | Protected`, so `Attributes[Catalan]` is `{Constant, Protected}`. Interned
name `SYM_Catalan` (`src/sym_names.c`), docstring in `src/info.c`. It carries no
DownValues, and `D[Catalan, x] -> 0` follows from `Constant`.

**Representation & numeric value.** A bare `EXPR_SYMBOL`, kept exact under evaluation.
`N[Catalan]` returns the tabled machine double; `N[Catalan, k]` runs
`fill_mpfr_catalan` → MPFR's `mpfr_const_catalan` at the requested precision, returning
an `EXPR_MPFR`. Machine value 0.915966;
`N[Catalan, 50] = 0.915965594177219015054603514932384110774149374281673`.

**Usage & limits.** `Catalan` arises as the closed form of various definite integrals
and sums; `NumericQ[Catalan]` is `True` (whitelisted in `is_numeric_quantity`,
`src/core.c`). Its irrationality is open; Mathilda treats it as an exact symbol whose
digits are produced on demand by `N`.
