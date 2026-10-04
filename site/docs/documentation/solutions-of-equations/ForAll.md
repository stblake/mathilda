# ForAll

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ForAll[x, expr]`**

The quantified statement that expr is True for all values of x. ForAll\[{x1, x2, ...}, expr\] binds several variables and ForAll\[x, cond, expr\] quantifies over values satisfying cond. ForAll is inert on its own (HoldAll); it is eliminated by Reduce or Resolve over the reals.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

Every real square is nonnegative, so True

```mathematica
In[1]:= Resolve[ForAll[x, x^2 >= 0]]
Out[1]= True
```

The parabola stays above zero exactly when a >= 0

```mathematica
In[2]:= Resolve[ForAll[x, a x^2 + 1 > 0], Reals]
Out[2]= a >= 0
```

Positive-definite on the discriminant band -2 < a < 2

```mathematica
In[3]:= Resolve[ForAll[x, x^2 + a x + 1 > 0], Reals]
Out[3]= -2 < a < 2
```

Below every square iff x <= 0

```mathematica
In[4]:= Reduce[ForAll[y, x <= y^2], x, Reals]
Out[4]= x <= 0
```

## Algorithm

reduce_qe.c

Quantifier elimination for `Reduce` (REDUCE_PLAN.md, Phase 7): the front-end

```text
for the `Exists`, `ForAll` and `Resolve` heads.  See reduce_qe.h for the shape
```

of the method and the three-case (by free-variable count) routing.

This file owns the front-end only -- quantifier normalisation (flatten a same-kind chain, fold a 3-argument condition), free-variable collection, the fully-quantified DECISION path (Case A, which reuses the whole Reduce engine), the routing to reduce_cad_qe for the parametric (>=1 free variable) path, and the recursive composition that eliminates an alternating quantifier prefix

```text
inner-block-first.  The CAD projection/lifting/fold and the multi-free-variable
```

emission live in reduce_cad.c.

Hard invariant: any decline (a malformed node, a non-Reals domain, or an undecidable/unsolvable sub-problem -- including an alternating sub-block whose inner elimination declines) returns NULL, leaving the input unevaluated -- never a wrong formula.

## Implementation notes

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

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Exists](../../solutions-of-equations/Exists/), [HoldAll](../../expression-information/HoldAll/), [Reduce](../../solutions-of-equations/Reduce/), [Resolve](../../solutions-of-equations/Resolve/)

- G. E. Collins and H. Hong, *Partial Cylindrical Algebraic Decomposition for Quantifier Elimination*, J. Symbolic Computation **12** (1991) 299-328.
- S. McCallum, *An Improved Projection Operation for Cylindrical Algebraic Decomposition*, in *Quantifier Elimination and Cylindrical Algebraic Decomposition* (Springer, 1998) 242-268.
- Source: [`src/solve/reduce_qe.c`](https://github.com/stblake/mathilda/blob/main/src/solve/reduce_qe.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)

## Notes & additional examples

### Notes

`ForAll[x, expr]` is the quantified statement that `expr` holds for *every* value
of `x`; `ForAll[{x1, x2, ...}, expr]` binds several variables and
`ForAll[x, cond, expr]` quantifies over values satisfying `cond` (read as
`cond` implies `expr`).

Like `Exists`, `ForAll` is **inert on its own** (`HoldAll`) and is evaluated only
when `Reduce` or `Resolve` eliminates it, over the **Reals**. A fully quantified
sentence returns `True` or `False`; a statement with free parameters returns the
condition on them — the classic use being the parameter band under which a
family of polynomials keeps one sign, as in the positive-definite quadratic
`x^2 + a x + 1 > 0` holding for all `x` exactly on `-2 < a < 2`.

Internally a universal is decided by the dual emptiness question to the
existential: `ForAll[x, g]` is true exactly when `!g` has no solution. Parametric
cases are projected with a Cylindrical Algebraic Decomposition, and alternating
prefixes are eliminated inner-block-first. An undecidable sign or an unsupported
construct leaves the statement unevaluated.
