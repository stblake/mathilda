# Xor

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Xor[e1, e2, ...]`**

The logical exclusive OR of the ei: True when an odd number of the arguments are True.  Folds literal Booleans and cancels duplicate arguments (a Xor a is False); Xor\[\] is False and Xor\[e\] is e.

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Flat`, `OneIdentity`, `Orderless`, `Protected`.

## References

**See also:** [Flat](../../expression-information/Flat/), [Orderless](../../expression-information/Orderless/), [OneIdentity](../../expression-information/OneIdentity/)

- Source: [`src/boolean.c`](https://github.com/stblake/mathilda/blob/main/src/boolean.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_boolean.c`](https://github.com/stblake/mathilda/blob/main/tests/test_boolean.c)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)
