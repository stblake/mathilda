---
source: src/core.c
---
**Algorithm.** `builtin_activate` (`src/core.c`) accepts exactly one argument (any other
arity returns `NULL`) and calls `activate_recursive`, which walks the expression tree
replacing every compound head `Inactive[h]` (an `Inactive` call of arity 1) by `h` and
rebuilding the node. The reactivated tree is handed back to the evaluator's fixed-point
loop, which then re-evaluates the now-active heads — so `Activate[Inactive[Integrate][g,
x]]` becomes `Integrate[g, x]` and actually integrates.

**Data structures.** Pure `Expr`-tree recursion. A non-function node is deep-copied
(`expr_copy`); a function node reactivates its head and every argument into a fresh
`Expr**` buffer that is passed to `expr_new_function` (which consumes the head and args).
The `Inactive[h] -> h` rewrite short-circuits before the generic head/argument walk, so
nested inert heads are peeled in one pass.

**Complexity / limits.** `O(size of expression)` for the single rewrite pass, plus
whatever re-evaluation the reactivated heads subsequently trigger. `Activate` is
`ATTR_PROTECTED` and is the exact inverse of `Inactive`.
