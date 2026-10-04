---
source: src/numeric.c
references:
  - "S. R. Finch, *Mathematical Constants* (Cambridge Univ. Press, 2003), §2.15."
---
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
