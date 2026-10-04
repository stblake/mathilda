# Exists

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Exists[x, expr]`**

The quantified statement that there exists a value of x for which expr is True.  Exists\[{x1, x2, ...}, expr\] binds several variables and Exists\[x, cond, expr\] restricts to values satisfying cond. Exists is inert on its own (HoldAll); it is eliminated by Reduce or Resolve over the reals.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

There is a real square root of 2, so the sentence is True

```mathematica
In[1]:= Resolve[Exists[x, x^2 == 2]]
Out[1]= True
```

No real x solves it, so False

```mathematica
In[2]:= Resolve[Exists[x, x^2 + 1 == 0]]
Out[2]= False
```

Solvable for x exactly when a != 0

```mathematica
In[3]:= Resolve[Exists[x, a x == 1], Reals]
Out[3]= a != 0
```

The image of squaring is the nonnegative half-line

```mathematica
In[4]:= Reduce[Exists[y, x == y^2], x]
Out[4]= x >= 0
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

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Reduce](../../solutions-of-equations/Reduce/), [Resolve](../../solutions-of-equations/Resolve/)

- G. E. Collins and H. Hong, *Partial Cylindrical Algebraic Decomposition for Quantifier Elimination*, J. Symbolic Computation **12** (1991) 299-328.
- S. McCallum, *An Improved Projection Operation for Cylindrical Algebraic Decomposition*, in *Quantifier Elimination and Cylindrical Algebraic Decomposition* (Springer, 1998) 242-268.
- Source: [`src/solve/reduce_qe.c`](https://github.com/stblake/mathilda/blob/main/src/solve/reduce_qe.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)

## Notes & additional examples

### Notes

`Exists[x, expr]` is the quantified statement that *some* value of `x` makes
`expr` true; `Exists[{x1, x2, ...}, expr]` binds several variables and
`Exists[x, cond, expr]` restricts the witness to values satisfying `cond`.

`Exists` is **inert on its own** — it carries `HoldAll`, keeps its bound
variables symbolic, and does not evaluate until `Reduce` or `Resolve` eliminates
it. Evaluate it with `Resolve[Exists[...]]`, or place it inside `Reduce[...]` to
get the condition on the remaining free variables. Quantified problems are solved
over the **Reals**; an explicit non-`Reals` domain declines.

When every variable is bound (a fully quantified *sentence*), the result is
`True` or `False`. When free variables remain — as in `Reduce[Exists[y, x ==
y^2], x]` — the bound variable `y` is projected away by a Cylindrical Algebraic
Decomposition and the answer is a condition on the free ones (here `x >= 0`).
Alternating quantifier prefixes such as `Exists[x, ForAll[y, ...]]` are
eliminated inner-block-first. An undecidable sign or an unsupported construct
leaves the statement unevaluated rather than guessed.
