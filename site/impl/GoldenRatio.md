---
source: src/numeric.c
references:
  - "S. R. Finch, *Mathematical Constants* (Cambridge Univ. Press, 2003), §1.2."
---
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
