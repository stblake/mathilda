---
references:
  - "G. E. Collins and H. Hong, *Partial Cylindrical Algebraic Decomposition for Quantifier Elimination*, J. Symbolic Computation **12** (1991) 299-328."
  - "S. McCallum, *An Improved Projection Operation for Cylindrical Algebraic Decomposition*, in *Quantifier Elimination and Cylindrical Algebraic Decomposition* (Springer, 1998) 242-268."
source: src/solve/reduce_qe.c
---
**Algorithm.** `Exists` is **inert on its own**: `builtin_exists` returns `NULL`
and the head carries `HoldAll`, so the bound variables stay symbolic and the node
simply holds a binding. The elimination happens when `Reduce` or `Resolve` meets
the head and calls `reduce_qe_dispatch`. That front-end runs `qe_normalize`, which
peels a maximal chain of same-kind quantifiers (collecting their bound-variable
names and conjoining any 3-argument restriction `Exists[x, cond, g]` into a side
condition), and then folds the restriction into the matrix with
`qe_apply_condition` — for `Exists` this is `cond && M`. A different-kind inner
quantifier makes the prefix *alternating*; it is eliminated inner-block-first by a
recursive `reduce_qe_dispatch` on the body, and only then is the outer `Exists`
block eliminated over the resulting quantifier-free formula.

With the bound block fixed, the number of **free** variables selects the method.
Zero free variables is a decision (`qe_decide`): `Exists[{B}, g]` is `True` unless
the solution set of `g` is empty, which it tests by re-entering `Reduce[g, B,
Reals]` and asking whether the result is `False` — the robust emptiness question,
so a tautological region that `Reduce` reports without literally simplifying to
`True` is still decided correctly. One or more free variables is parametric QE
(`qe_parametric`): after the same `Reals` piecewise/radical preprocessing the base
engine uses, the matrix becomes an `RForm` handed to `reduce_cad_qe`, which runs an
iterated McCallum projection over all variables and emits a sign formula over the
free-variable subspace only.

**Data structures.** `Expr` trees throughout; bound and free variable names are
interned-pointer arrays; the matrix is an `RForm` (DNF of polynomial relations).
The CAD projection stacks (`PolySet`) and cell tree (`CADRegion`) live in
`reduce_cad.c` and use the exact `qqbar` real-algebraic sign oracle.

**Complexity / limits.** Quantified problems are supported over the **Reals**
only; an explicit non-`Reals` domain declines. Cost is that of the CAD — doubly
exponential in the variable count. The hard invariant is soundness: any malformed
node, non-`Reals` domain, or undecidable / unsupported sub-problem (including an
alternating inner block whose elimination declines) returns `NULL`, so the
`Exists` statement is left unevaluated rather than answered wrongly.
