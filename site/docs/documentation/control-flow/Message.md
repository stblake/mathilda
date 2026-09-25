# Message

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Message[sym::tag, args] prints the named message unless messages are suppressed.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

- `HoldFirst`.
- Returns `Null`.

**Attributes:** `HoldFirst`, `Protected`.

## References

**See also:** [Check](../../control-flow/Check/), [HoldFirst](../../other-advanced/HoldFirst/)

- Source: [`src/message.c`](https://github.com/stblake/mathilda/blob/main/src/message.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_parallelmixedtower.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parallelmixedtower.c)
