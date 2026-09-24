# Check

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Check[expr, failexpr] returns failexpr if a message is generated while evaluating expr, otherwise the value of expr.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

- `HoldAll`. A `Throw` inside `expr` propagates (it is not a message).
- Typically wrapped as `Quiet[Check[expr, failexpr]]` to detect a failure without printing its message.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Throw](../../control-flow/Throw/)

- Source: [`src/message.c`](https://github.com/stblake/mathilda/blob/main/src/message.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_parallelmixedtower.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parallelmixedtower.c)
