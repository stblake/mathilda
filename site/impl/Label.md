---
source: src/funcprog.c
---
**Algorithm.** `Label[tag]` marks a jump target inside a `CompoundExpression`.
`Label` is `Protected`. `builtin_label` validates arity (exactly one argument) and
returns `expr_new_symbol(SYM_Null)`, so evaluated as an ordinary statement a
`Label` is a no-op worth `Null`. Its real role is passive: the *raw held*
`Label[tag]` node, as it appears literally in the enclosing `CompoundExpression`,
is what `builtin_compoundexpression` scans for when it consumes a `Goto[tag]`
sentinel (see [Goto](Goto.md)). Tags are compared structurally (conventionally a
literal symbol or integer).

**Data structures.** None of its own; it is a plain two-node `EXPR_FUNCTION` that
`CompoundExpression` reads by position among its statements.

**Complexity / limits.** A `Label` is meaningful only as an explicit element of a
`CompoundExpression`; the matching `Goto` resolves it by a linear scan of that
compound expression's statements, then of enclosing ones. A `Label` reached in
normal top-to-bottom flow simply evaluates to `Null` and execution continues.
