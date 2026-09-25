# Quiet

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Quiet[expr] evaluates expr with messages suppressed and returns its value. Quiet[expr, spec] suppresses all messages regardless of spec.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

- `HoldAll`, so the argument is evaluated under the suppression, not before it.
- A message still *fires* while suppressed (an enclosing `Check` sees it); only the printing is silenced.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Check](../../control-flow/Check/)

- Source: [`src/message.c`](https://github.com/stblake/mathilda/blob/main/src/message.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_parallelmixedtower.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parallelmixedtower.c)
