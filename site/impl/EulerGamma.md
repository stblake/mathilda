---
source: src/numeric.c
references:
  - "S. R. Finch, *Mathematical Constants* (Cambridge Univ. Press, 2003), §1.5."
---
**Definition.** `EulerGamma` is the Euler–Mascheroni constant γ ≈ 0.5772156649, the
limit `lim_{n->Infinity} (HarmonicNumber[n] - Log[n])` and the constant term of the
Laurent expansion of `Zeta[s]` about `s = 1`. It is a constant *symbol*: row
`{ "EulerGamma", 0.5772156649…, fill_mpfr_eulergamma }` of `kConstants[]` in
`src/numeric.c` supplies its value and the `numericalize_init` loop stamps it
`Constant | Protected`; the same attributes are *additionally* stamped by
`eulergamma_init` in `src/special_functions/eulergamma.c`, which owns nothing but the
symbol's identity. So `Attributes[EulerGamma]` is `{Constant, Protected}`. Interned
name `SYM_EulerGamma` (`src/sym_names.c`), docstring in `src/info.c`. `EulerGamma`
has no DownValues; it surfaces as the constant term of `Series[Zeta[s], {s, 1, n}]`,
and `D[EulerGamma, x] -> 0` follows from `Constant`.

**Representation & numeric value.** A bare `EXPR_SYMBOL`, kept exact through symbolic
computation. `N[EulerGamma]` returns the tabled machine double; `N[EulerGamma, k]`
runs `fill_mpfr_eulergamma` → MPFR's `mpfr_const_euler` at the requested precision,
returning an `EXPR_MPFR`. Machine value 0.577216;
`N[EulerGamma, 50] = 0.577215664901532860606512090082402431042159335939923`.

**Usage & limits.** `EulerGamma` appears in the Zeta Laurent expansion, PolyGamma /
digamma values and the Stieltjes-constant family; `NumericQ[EulerGamma]` is `True`
(whitelisted in `is_numeric_quantity`, `src/core.c`). Whether γ is irrational is
famously open; Mathilda simply treats it as an exact symbol whose digits are produced
on demand by `N`.
