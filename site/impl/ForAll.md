---
references:
  - "G. E. Collins and H. Hong, *Partial Cylindrical Algebraic Decomposition for Quantifier Elimination*, J. Symbolic Computation **12** (1991) 299-328."
  - "S. McCallum, *An Improved Projection Operation for Cylindrical Algebraic Decomposition*, in *Quantifier Elimination and Cylindrical Algebraic Decomposition* (Springer, 1998) 242-268."
source: src/solve/reduce_qe.c
---
**Algorithm.** `ForAll` is **inert on its own**: `builtin_forall` returns `NULL`
and the head carries `HoldAll`, so the bound variables stay symbolic until `Reduce`
or `Resolve` meets the head and calls `reduce_qe_dispatch`. `qe_normalize` peels a
maximal chain of same-kind `ForAll`s, collects their bound variables, and conjoins
any 3-argument restriction `ForAll[x, cond, g]` into a side condition.
`qe_apply_condition` then folds the restriction as a bounded universal:
`ForAll[x, cond, M]` becomes `!cond || M`. An inner different-kind quantifier makes
the prefix alternating and is eliminated inner-block-first by a recursive dispatch
before the outer `ForAll` block is eliminated over the resulting quantifier-free
formula.

The free-variable count then selects the method. Zero free variables is a decision
(`qe_decide`): `ForAll[{B}, g]` is `True` exactly when `!g` is unsatisfiable, which
it tests by re-entering `Reduce[!g, B, Reals]` and checking that the result is
`False` — so the universal is decided through the same robust emptiness query as
the existential, just on the negated matrix. One or more free variables is
parametric QE (`qe_parametric` → `reduce_cad_qe`): the matrix is preprocessed for
`Reals` (selector case-splitting, radical rationalisation), turned into an `RForm`,
and run through an iterated McCallum projection whose sign formula is emitted over
the free-variable subspace — producing, e.g., the discriminant band under which a
parametric quadratic is positive for all `x`.

**Data structures.** `Expr` trees; interned-pointer name arrays for bound/free
variables; the matrix as an `RForm` (DNF of polynomial relations with `R_EQ` /
`R_LT` / `R_LE` codes). CAD projection stacks (`PolySet`), the `CADRegion` cell
tree, and the exact `qqbar` sign oracle live in `reduce_cad.c`.

**Complexity / limits.** Supported over the **Reals** only — a non-`Reals` domain
declines — with CAD cost doubly exponential in the variable count. Soundness is the
hard invariant: a malformed node, a non-`Reals` domain, an undecidable sign, or an
unsupported / declining sub-block returns `NULL`, leaving the `ForAll` statement
unevaluated instead of guessed.
