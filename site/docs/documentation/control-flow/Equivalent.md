# Equivalent

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Equivalent[e1, e2, ...]`**

The logical equivalence e1 \\[Equivalent\] e2 \\[Equivalent\] ...: True when all of the ei have the same truth value.  Folds literal Booleans and cancels duplicate arguments; Equivalent\[\] and Equivalent\[e\] are True.

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Flat`, `OneIdentity`, `Orderless`, `Protected`.

## References

**See also:** [Flat](../../expression-information/Flat/), [Orderless](../../expression-information/Orderless/), [OneIdentity](../../expression-information/OneIdentity/), [LogicalExpand](../../solutions-of-equations/LogicalExpand/), [Reduce](../../solutions-of-equations/Reduce/), [FindInstance](../../solutions-of-equations/FindInstance/)

- Source: [`src/boolean.c`](https://github.com/stblake/mathilda/blob/main/src/boolean.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_boolean.c`](https://github.com/stblake/mathilda/blob/main/tests/test_boolean.c)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)
