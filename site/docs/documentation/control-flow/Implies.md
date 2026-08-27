# Implies

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Implies[p, q]`**

The material implication p \\[Implies\] q, equivalent to !p || q. Implies\[False, q\] and Implies\[p, True\] are True, Implies\[True, q\] is q, and Implies\[p, False\] is !p.

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [LogicalExpand](../../solutions-of-equations/LogicalExpand/), [Reduce](../../solutions-of-equations/Reduce/)

- Source: [`src/boolean.c`](https://github.com/stblake/mathilda/blob/main/src/boolean.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_boolean.c`](https://github.com/stblake/mathilda/blob/main/tests/test_boolean.c)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)
