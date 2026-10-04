---
source: src/numeric.c
---
**Definition.** `E` is Euler's number *e* = 2.71828…, the base of the natural
logarithm and the unique real with `D[E^x, x] == E^x`. It is a constant *symbol*:
row `{ "E", M_E, fill_mpfr_e }` of the `kConstants[]` table in `src/numeric.c` supplies
its numeric value, and the `numericalize_init` loop stamps it `Constant | Protected`
(so `Attributes[E]` is `{Constant, Protected}`; the base `Protected` flag is also set
in the default-attributes table in `src/attr.c`). Interned name `SYM_E`
(`src/sym_names.c`), docstring in `src/info.c`. `E` has no DownValues; the identities
that fold back to it — `Log[E] -> 1`, `Log[E^3] -> 3`, `Exp[1] -> E`, and the series /
limit characterisations `Sum[1/n!, {n,0,Infinity}] -> E`,
`Limit[(1+1/n)^n, n->Infinity] -> E` — are produced by the log-exp, sum and limit
subsystems. `D[E, x] -> 0` follows from `Constant`.

**Representation & numeric value.** `E` is a bare `EXPR_SYMBOL` that prints as `E` and
stays exact under evaluation. `N[E]` returns the machine double `M_E`; `N[E, k]` runs
`fill_mpfr_e`, which evaluates `mpfr_exp(1)` at the requested precision (MPFR ships no
dedicated `const_e`) and returns an `EXPR_MPFR`. Machine value 2.71828;
`N[E, 50] = 2.71828182845904523536028747135266249775724709369996`.

**Usage & limits.** `E` is the fixed point of the exponential closed forms and is
recognised as numeric (`NumericQ[E]` is `True`, via `is_numeric_quantity` in
`src/core.c`). Being transcendental it is never rationalised; its digits are produced
on demand by `N`.
