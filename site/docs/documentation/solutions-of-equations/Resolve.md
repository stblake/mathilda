# Resolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Resolve[expr]`**

**`Resolve[expr, dom]`**

Eliminates the quantifiers (Exists, ForAll) from expr over the domain dom (Reals; the default and only supported domain), returning an equivalent quantifier-free statement -- True or False for a fully quantified sentence, or a condition on the remaining free variables. Parametric elimination is supported for a single free variable; an undecidable, alternating, or higher-dimensional case is left unevaluated rather than guessed.

## Examples

_No verified examples yet for this function._

## Algorithm

reduce_qe.c

Quantifier elimination for `Reduce` (REDUCE_PLAN.md, Phase 7): the front-end

```text
for the `Exists`, `ForAll` and `Resolve` heads.  See reduce_qe.h for the shape
```

of the method and the three-case (by free-variable count) routing.

This file owns the front-end only -- quantifier normalisation (flatten a same-kind chain, fold a 3-argument condition), free-variable collection, the fully-quantified DECISION path (Case A, which reuses the whole Reduce engine), and the routing to reduce_cad_qe for the parametric single-free-variable path

```text
(Case B).  The CAD projection/lifting/fold machinery lives in reduce_cad.c.
```

Hard invariant: any decline (a malformed node, an alternating quantifier prefix, >=2 free variables, a non-Reals domain, or an undecidable/unsolvable sub-problem) returns NULL, leaving the input unevaluated -- never a wrong formula.

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Exists](../../solutions-of-equations/Exists/), [ForAll](../../solutions-of-equations/ForAll/)

- Source: [`src/solve/reduce_qe.c`](https://github.com/stblake/mathilda/blob/main/src/solve/reduce_qe.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)
