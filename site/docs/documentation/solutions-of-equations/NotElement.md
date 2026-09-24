# NotElement

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NotElement[x, dom]`**

The statement that x is not an element of the domain dom -- the negation of Element\[x, dom\].  Decides to True or False when the membership decides, and stays symbolic otherwise.

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [LogicalExpand](../../solutions-of-equations/LogicalExpand/)

- Source: [`src/solve/reduce_companions.c`](https://github.com/stblake/mathilda/blob/main/src/solve/reduce_companions.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_deriv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv.c)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)
