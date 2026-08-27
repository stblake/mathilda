# NotElement

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NotElement[x, dom]`**

The statement that x is not an element of the domain dom -- the negation of Element\[x, dom\].  Decides to True or False when the membership decides, and stays symbolic otherwise.

## Examples

_No verified examples yet for this function._

## Algorithm

reduce_companions.c

```text
Companion builtins for `Reduce` (REDUCE_PLAN.md, Phase 8).  v1: LogicalExpand
+ a minimal NotElement head.
```

LogicalExpand distributes a logical statement to disjunctive normal form (an Or of Ands of literals), applying idempotence / complementation / absorption contractions, and collapsing to True (tautology) or False (contradiction)

```text
when the statement decides.  Every non-connective subexpression is treated as
```

an OPAQUE Boolean atom -- a symbol, a relation `x == a`, a membership

```text
`Element[..]` -- with NO domain reasoning, exactly as Mathematica's
LogicalExpand does.  Two relational atoms are complementary iff one is the
```

(head-flipped) logical negation of the other (`x==a` / `x!=a`, `x<1` / `x>=1`,

```text
`Element` / `NotElement`, `a` / `!a`).
```

The True/False collapse is sound *and* complete without truth-table enumeration: over independent opaque atoms a DNF is unsatisfiable iff every clause holds a complementary pair -- i.e. it distributes to ZERO surviving

```text
clauses.  So `phi` empty => False, and `Not[phi]` empty => True (the negation
```

is unsatisfiable, hence `phi` is a tautology).

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [LogicalExpand](../../solutions-of-equations/LogicalExpand/)

- Source: [`src/solve/reduce_companions.c`](https://github.com/stblake/mathilda/blob/main/src/solve/reduce_companions.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)
