---
source: src/calculus/deriv.c
---
**Algorithm.** `Dt` shares the native differentiation core with `D`
(`compute_deriv` in src/calculus/deriv.c). `builtin_dt` handles two modes.
`Dt[f]` (one argument) computes the **total derivative**: it calls
`compute_deriv(f, NULL, NULL)` with a NULL differentiation variable, so unknown
symbols are *not* treated as constants — each contributes a `Dt[sym]`
differential term, and the usual product/quotient/chain rules thread through.
`Dt[f, var, ...]` is defined to be identical to `D[f, var, ...]` (the partial
derivative) and is forwarded to the same per-spec loop used by `builtin_d`
(`parse_var_spec` + `higher_order_partial` / `array_higher_order` /
`compute_deriv_symbolic_order`). Malformed specs emit a `D::dvar`-style message
and return unevaluated.

**Data structures.** `Expr*` tree transformation only; results are returned
un-reduced and folded by the outer evaluator.

**Complexity / limits.** Linear per pass in the tree size. The total-derivative
mode distinguishes itself from `D` solely by the NULL variable that disables the
constant short-circuit; everything else (rules, ownership, fixed-point folding)
is shared with `D`.
