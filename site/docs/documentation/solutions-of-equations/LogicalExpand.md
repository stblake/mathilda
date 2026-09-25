# LogicalExpand

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LogicalExpand[expr]`**

Expands the logical combination expr -- of equations, inequalities and Boolean atoms -- into disjunctive normal form (an Or of Ands), applying distributive, De Morgan, idempotence, complementation and absorption laws, and expanding Implies and Xor.  Returns True for a tautology and False for a contradiction.  Every non-logical subexpression is treated as an opaque Boolean atom (no domain reasoning), so e.g. x==a and x!=a are complementary literals.

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Implies](../../control-flow/Implies/), [Xor](../../control-flow/Xor/), [Reduce](../../solutions-of-equations/Reduce/), [Element](../../simplification/Element/), [NotElement](../../solutions-of-equations/NotElement/)

- Source: [`src/solve/reduce_companions.c`](https://github.com/stblake/mathilda/blob/main/src/solve/reduce_companions.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)
