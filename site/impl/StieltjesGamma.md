---
source: src/special_functions/stieltjesgamma.c
references:
  - "DLMF §25.2.4 — the Stieltjes constants in the Laurent expansion of zeta about s = 1."
---
**Algorithm.** `builtin_stieltjesgamma` represents the Stieltjes constants
`StieltjesGamma[n] = gamma_n`, the coefficients of the Laurent expansion of the
Riemann zeta function about `s = 1`:
`zeta(s) = 1/(s-1) + Sum_{n>=0} ((-1)^n/n!) gamma_n (s-1)^n`. It is deliberately
**inert**: the only reduction is `StieltjesGamma[0] -> EulerGamma`; every other
order (and all generic symbolic behaviour) is left to the evaluator. It is the
natural output of `Series[Zeta[x], {x, 1, n}]`. The module owns only the symbol's
identity and the `n = 0` reduction; the docstring lives in `info.c`.

**Data structures.** `Expr` only — there is no numeric backend. No ND kernel and
not on `pack.c`'s `AWARE` list; `Compile[]` does not lower it
(`CompileDiagnostics` reports `Compiled -> False`). Attributes: `Listable`,
`Protected` (notably **not** `NumericFunction` — the higher constants have no
numeric evaluation here).

**Complexity / limits.** `O(1)` symbolic reduction. The higher Stieltjes
constants (`n >= 1`) have no elementary closed form and are not evaluated
numerically; they stay symbolic.
