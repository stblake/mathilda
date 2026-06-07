---
source: src/calculus/deriv.c
---
**Algorithm.** `D[f, x]` is computed by a native, dispatch-driven C
differentiator (`builtin_d` → `compute_deriv` in src/calculus/deriv.c), which
replaced the old rule-based `src/internal/deriv.m` (now a no-op stub kept only
so users can drop in custom `Dt` identities). `builtin_d` first splits trailing
arguments into `NonConstants -> ...` options and variable specs, then applies
each spec sequentially so that mixed partials `D[f, x, y]` and higher orders
`D[f, {x, n}]` and array forms `D[f, {{x1,...,xN}}]` all reduce to repeated
single-variable differentiation (`higher_order_partial`, `array_higher_order`,
and `compute_deriv_symbolic_order` for symbolic order `n`).

The core `compute_deriv(f, x, nonconsts)` performs a single head-symbol dispatch
per node: linearity for `Plus`; the general n-factor **product rule** for
`Times`; `Power` handled with the combined power/exponential/general
`u^v` rule; the **quotient** case falls out of `Power[..,-1]`; and elementary
unary functions (Sin, Cos, Exp, Log, ArcTan, …) get their derivative from a
table lookup (`elementary_fprime`) composed with the **chain rule**. Unknown
single- and multi-argument functions get the generic chain rule
(`chain_rule_unknown`), emitting `Derivative[...][f]` factors. Constant
subtrees are short-circuited by a tailored structural walk (`expr_free_of`)
rather than calling the generic `FreeQ` builtin — this is the main speedup over
the rule-based version. Each builder returns plain un-reduced trees (e.g.
`Plus[0, x]`, `Times[1, x]`); the outer Mathilda fixed-point evaluator folds the
arithmetic.

**Data structures.** Pure `Expr*` tree transformation; no auxiliary numeric
representation. Symbolic-order derivatives produce closed forms where possible
and otherwise fall back to an unevaluated `D[...]`.

**Complexity / limits.** Linear in the size of the input tree per
differentiation pass (with constant-subtree pruning); `D[f, x, n]` costs n
passes. `NonConstants` is honoured for ordinary specs but not threaded through
the closed-form symbolic-order path.
