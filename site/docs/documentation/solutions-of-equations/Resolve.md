# Resolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Resolve[expr]`**

**`Resolve[expr, dom]`**

Eliminates the quantifiers (Exists, ForAll) from expr over the domain dom (Reals; the default and only supported domain), returning an equivalent quantifier-free statement -- True or False for a fully quantified sentence, or a condition on the remaining free variables (one or more). Alternating quantifier prefixes are eliminated inner-block-first. An undecidable sign, a non-Reals domain, or an unsupported case is left unevaluated rather than guessed.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

An inverse for a exists exactly when a != 0

```mathematica
In[1]:= Resolve[Exists[x, a x == 1], Reals]
Out[1]= a != 0
```

Discriminant band for a positive-definite quadratic

```mathematica
In[2]:= Resolve[ForAll[x, x^2 + a x + 1 > 0], Reals]
Out[2]= -2 < a < 2
```

Fully quantified sentence, decides to True

```mathematica
In[3]:= Resolve[Exists[x, x^2 == 2 && x > 0]]
Out[3]= True
```

An alternating prefix, eliminated inner-block-first

```mathematica
In[4]:= Resolve[ForAll[x, Exists[y, y > x]]]
Out[4]= True
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

**Algorithm.** `Resolve` is the dedicated quantifier-elimination front-end.
`builtin_resolve` accepts one or two arguments, requires the first to be a
quantified head (`Exists` / `ForAll`) — otherwise it declines — and calls
`reduce_qe_dispatch(qexpr, dom)`, the same engine `Reduce` invokes when its input
happens to be top-level quantified. The dispatch normalises the quantifier prefix
(`qe_normalize` peels a maximal same-kind chain, records the bound variables, and
carries any 3-argument restriction as a side condition folded in by
`qe_apply_condition`: `cond && M` for `Exists`, `!cond || M` for `ForAll`). An
alternating prefix is eliminated **inner-block-first** by recursion: the inner
different-kind quantifier is resolved to a quantifier-free `psi`, this block's
restriction is folded onto it, and the outer block is then re-eliminated over the
non-alternating result — composing to arbitrary alternation depth.

Elimination of a single block branches on the number of remaining free variables.
Zero free variables is a decision via `qe_decide`, which re-enters `Reduce` on the
matrix (or its negation, for `ForAll`) and reads emptiness, returning `True` or
`False`. One or more free variables is parametric QE via `qe_parametric` →
`reduce_cad_qe`: after `Reals` preprocessing the matrix becomes an `RForm` fed to an
iterated McCallum projection, whose sign formula is emitted over the free-variable
subspace (free-variable names are sorted alphabetically first, so the answer is
stable). A fully quantified sentence therefore resolves to `True`/`False`; a
partially quantified one resolves to a condition on its free variables.

**Data structures.** `Expr` trees; interned-pointer arrays for the bound and free
variable names; the matrix as an `RForm`. The CAD machinery (`PolySet` projection
stacks, the `CADRegion` cell tree, the `qqbar` sign oracle) lives in
`reduce_cad.c`.

**Complexity / limits.** `Reals` is the default and only supported domain — any
other declines. Cost is that of the CAD (doubly exponential in the variable count).
By the engine's hard invariant, an undecidable sign, a non-`Reals` domain, an
unsupported construct, or an inner block whose elimination declines all return
`NULL`, leaving the `Resolve[...]` call unevaluated rather than returning a wrong
formula.

**Attributes:** `Protected`.

## References

**See also:** [Exists](../../solutions-of-equations/Exists/), [ForAll](../../solutions-of-equations/ForAll/), [Abs](../../arithmetic/Abs/), [Min](../../data-structures/Min/), [Max](../../data-structures/Max/), [Piecewise](../../control-flow/Piecewise/), [Sign](../../arithmetic/Sign/), [UnitStep](../../elementary-functions/UnitStep/)

- G. E. Collins and H. Hong, *Partial Cylindrical Algebraic Decomposition for Quantifier Elimination*, J. Symbolic Computation **12** (1991) 299-328.
- S. McCallum, *An Improved Projection Operation for Cylindrical Algebraic Decomposition*, in *Quantifier Elimination and Cylindrical Algebraic Decomposition* (Springer, 1998) 242-268.
- Source: [`src/solve/reduce_qe.c`](https://github.com/stblake/mathilda/blob/main/src/solve/reduce_qe.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)

## Notes & additional examples

### Notes

`Resolve[expr]` eliminates the quantifiers (`Exists`, `ForAll`) from `expr`,
returning an equivalent quantifier-free statement: `True` or `False` for a fully
quantified sentence, or a condition on the remaining free variables otherwise. It
is the dedicated quantifier-elimination front-end — its argument must be a
quantified head, and it shares the engine `Reduce` uses when its own input happens
to be quantified. `Resolve[expr, dom]` names the domain, but **Reals** is the
default and only supported one.

Alternating quantifier prefixes such as `ForAll[x, Exists[y, y > x]]` are handled
by eliminating the innermost block to a quantifier-free formula and then
re-eliminating the enclosing block over it, composing to arbitrary depth. The
parametric case is driven by a Cylindrical Algebraic Decomposition that projects
away the bound variables and emits a sign formula over the free ones.

`Resolve` follows the engine's soundness rule: an undecidable sign, a non-`Reals`
domain, or a construct outside the supported fragment leaves the call unevaluated
rather than returning a wrong formula.
