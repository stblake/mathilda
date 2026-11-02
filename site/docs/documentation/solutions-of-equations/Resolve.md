# Resolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Resolve[expr]`**

**`Resolve[expr, dom]`**

Eliminates the quantifiers (Exists, ForAll) from expr over the domain dom (Reals; the default and only supported domain), returning an equivalent quantifier-free statement -- True or False for a fully quantified sentence, or a condition on the remaining free variables (one or more). Alternating quantifier prefixes are eliminated inner-block-first. An undecidable sign, a non-Reals domain, or an unsupported case is left unevaluated rather than guessed.

## Examples (16)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= Resolve[ForAll[x, x^2 + a x + 1 > 0], Reals]
Out[1]= -2 < a < 2

In[2]:= Resolve[Exists[x, x^2 + b x + c == 0], Reals]
Out[2]= c <= 1/4 b^2

In[3]:= Resolve[ForAll[x, Abs[x] < d, x^2 < 9], Reals]
Out[3]= d <= 3

In[4]:= Resolve[ForAll[eps, eps > 0, Exists[del, del > 0, ForAll[x, 0 < Abs[x - 2] < del, Abs[(3 x - 1) - 5] < eps]]], Reals]
Out[4]= True

In[5]:= Resolve[ForAll[eps, eps > 0, Exists[nn, nn > 0, ForAll[x, x > nn, Abs[1/x] < eps]]], Reals]
Out[5]= True

In[6]:= Resolve[Exists[{x, y, z}, x + y + z == 1 && x >= 0 && y >= 0 && z >= 0 && x^2 + y^2 + z^2 > 4], Reals]
Out[6]= False
```

### Applications (10)

An inverse for a exists exactly when a != 0

```mathematica
In[7]:= Resolve[Exists[x, a x == 1], Reals]
Out[7]= a != 0
```

Discriminant band for a positive-definite quadratic

```mathematica
In[8]:= Resolve[ForAll[x, x^2 + a x + 1 > 0], Reals]
Out[8]= -2 < a < 2
```

Parametric elimination of a quartic: q below the minimum value 3 (p/4)^(4/3)

```mathematica
In[9]:= Resolve[Exists[x, x^4 + p x + q < 0], Reals]
Out[9]= p < 0 && q < (27/256 p^4)^(1/3) || p == 0 && q < 0 || p > 0 && q < (27/256 p^4)^(1/3)
```

An alternating prefix, eliminated inner block first

```mathematica
In[10]:= Resolve[ForAll[x, Exists[y, y > x]]]
Out[10]= True
```

A restricted quantifier whose Abs is case-split before the CAD

```mathematica
In[11]:= Resolve[ForAll[x, Abs[x] < d, x^2 < 9], Reals]
Out[11]= d <= 3
```

The epsilon-delta proof that 3 x - 1 tends to 5 at x = 2

```mathematica
In[12]:= Resolve[ForAll[eps, eps > 0, Exists[del, del > 0, ForAll[x, 0 < Abs[x - 2] < del, Abs[(3 x - 1) - 5] < eps]]], Reals]
Out[12]= True
```

The wrong limit 6 is refuted

```mathematica
In[13]:= Resolve[ForAll[eps, eps > 0, Exists[del, del > 0, ForAll[x, 0 < Abs[x - 2] < del, Abs[(3 x - 1) - 6] < eps]]], Reals]
Out[13]= False
```

A limit at infinity, decidable once 1/x is cleared to polynomial form

```mathematica
In[14]:= Resolve[ForAll[eps, eps > 0, Exists[nn, nn > 0, ForAll[x, x > nn, Abs[1/x] < eps]]], Reals]
Out[14]= True
```

A three-variable sentence decided False

```mathematica
In[15]:= Resolve[Exists[{x, y, z}, x + y + z == 1 && x >= 0 && y >= 0 && z >= 0 && x^2 + y^2 + z^2 > 4], Reals]
Out[15]= False
```

Transcendental body: no polynomial CAD, so it is left unevaluated

```mathematica
In[16]:= Resolve[Exists[x, Sin[x] == a], Reals]
Out[16]= Resolve[Exists[x, Sin[x] == a], Reals]
```

## Algorithm

reduce_qe.c

Quantifier elimination for `Reduce` (docs/notes/REDUCE_PLAN.md, Phase 7): the front-end

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
`reduce_cad_qe`. Before the CAD sees the matrix, `qe_parametric` runs
`reduce_piecewise_preprocess`, the real-domain preprocessing of the base `Reduce`
engine. It case-splits `Abs`, `Max`/`Min`, `Piecewise`, `Sign` and the other
selectors, rationalises square roots, and (since v0.289) clears variable
denominators to sign-exact polynomial atoms. A chained elimination can therefore carry
an intermediate rational atom such as `m <= 1/del` into the next level, which is what
lets ε–δ limit statements, limits at infinity and infinite limits decide. The
matrix then becomes an `RForm` fed to an iterated McCallum projection with the
bound variables projected out first. Its sign formula is emitted over the
free-variable subspace, with free-variable names sorted alphabetically so the answer
is stable. The CAD is built through `cad_build_augmented` (v0.321), the same
well-orientedness augment-and-retry loop as plain `Reduce`: a positive-dimensional
fibre nullification adds coefficients and rebuilds, instead of declining
immediately. A fully quantified sentence therefore resolves to `True`/`False`; a
partially quantified one resolves to a condition on its free variables.

**Data structures.** `Expr` trees; interned-pointer arrays for the bound and free
variable names; the matrix as an `RForm`. The CAD machinery (`PolySet` projection
stacks, the `CADRegion` cell tree, the `qqbar` sign oracle) lives in
`reduce_cad.c`.

**Complexity / limits.** `Reals` is the default and only supported domain — any
other declines, and a first argument that is not quantified declines too
(`Resolve[x > 0]` stays unevaluated). Cost is that of the CAD, doubly
exponential in the variable count. The CAD's loop-head deadline polls make
`TimeConstrained[Resolve[...], t]` return promptly. A transcendental body
(`Exists[x, Sin[x] == a]`) has no polynomial CAD and declines.
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

- `Resolve[expr]` and `Resolve[expr, Reals]` eliminate the quantifiers
  `Exists` and `ForAll` from `expr`. A fully quantified sentence gives
  `True` or `False`; otherwise the result is an equivalent quantifier-free
  condition on the free variables. The first argument must be quantified,
  so `Resolve[x > 0]` stays unevaluated. A quantified statement inside
  `Reduce` goes through the same engine.
- **Domain.** `Reals` is the default and the only domain supported. Any
  other domain declines.
- **What it handles.** It handles polynomial bodies, and also bodies with
  `Abs`, `Max`/`Min`, piecewise selectors, square roots, or rational atoms
  with variable denominators, all rewritten to polynomial branches first.
  The restricted form `ForAll[x, cond, body]` can appear at any level of
  an alternating prefix. Together these make ε–δ statements decidable for
  semialgebraic functions, including limits at infinity and infinite
  limits.
- **Soundness rule.** An undecidable sign, a non-`Reals` domain, a
  transcendental body, or a free leading coefficient that mixes a degree
  drop with an irrational boundary all leave the call unevaluated, rather
  than return a formula that might be wrong.
- `Protected`. Symbolic only.
- See also: `Reduce`, `Exists`, `ForAll`, `FindInstance`,
  `CylindricalDecomposition`.
