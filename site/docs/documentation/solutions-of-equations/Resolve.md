# Resolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Resolve[expr]`**

**`Resolve[expr, dom]`**

Eliminates the quantifiers (Exists, ForAll) from expr over the domain dom (Reals; the default and only supported domain), returning an equivalent quantifier-free statement -- True or False for a fully quantified sentence, or a condition on the remaining free variables (one or more). Alternating quantifier prefixes are eliminated inner-block-first. An undecidable sign, a non-Reals domain, or an unsupported case is left unevaluated rather than guessed.

## Examples

_No verified examples yet for this function._

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

**Attributes:** `Protected`.

## References

**See also:** [Exists](../../solutions-of-equations/Exists/), [ForAll](../../solutions-of-equations/ForAll/)

- Source: [`src/solve/reduce_qe.c`](https://github.com/stblake/mathilda/blob/main/src/solve/reduce_qe.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)
