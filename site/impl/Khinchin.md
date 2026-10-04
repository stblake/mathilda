---
source: src/numeric.c
references:
  - "D. H. Bailey, J. M. Borwein and R. E. Crandall, *On the Khintchine constant*, Math. Comp. **66** (1997) 417-431."
  - "S. R. Finch, *Mathematical Constants* (Cambridge Univ. Press, 2003), §1.8."
---
**Definition.** `Khinchin` is Khinchin's constant *K* (also "Khintchine's constant")
≈ 2.6854520011: the geometric mean of the partial quotients in the continued-fraction
expansion of almost every real number, `Product_{s>=1} (1 + 1/(s (s + 2)))^Log2[s]`.
It is a constant *symbol*: row `{ "Khinchin", 2.6854520011…, fill_mpfr_khinchin }` of
`kConstants[]` in `src/numeric.c` supplies its value and the `numericalize_init` loop
stamps it `Constant | Protected`, so `Attributes[Khinchin]` is `{Constant, Protected}`.
Interned name `SYM_Khinchin` (`src/sym_names.c`), docstring in `src/info.c`. No
DownValues; `D[Khinchin, x] -> 0` follows from `Constant`.

**Representation & numeric value.** A bare `EXPR_SYMBOL`, kept exact under evaluation.
`N[Khinchin]` returns the tabled machine double. `N[Khinchin, k]` runs
`fill_mpfr_khinchin`, which uses the geometrically convergent **Bailey–Borwein–Crandall
zeta series** `Log[K] Log[2] = Sum_{n>=1} ((Zeta[2n] - 1)/n) Sum_{k=1}^{2n-1}
(-1)^(k+1)/k`: `Zeta[2n]-1 ~ 4^-n`, so ~`bits/2` terms (each via `mpfr_zeta_ui`, with a
running alternating-harmonic inner sum) reach the target; `K = Exp[S / Log[2]]` at
`bits + 64` guard precision. Machine value 2.68545;
`N[Khinchin, 50] = 2.68545200106530644530971483548179569382038229399446`.

**Usage & limits.** `Khinchin` is the metric constant of continued fractions;
`NumericQ[Khinchin]` is `True` (whitelisted in `is_numeric_quantity`, `src/core.c`).
Whether *K* is irrational is open. The implementation limit mirrors `Glaisher`'s: the
filler is a bespoke convergent series rather than a library constant, with a term
cutoff at `2^-(bits+32)`.
