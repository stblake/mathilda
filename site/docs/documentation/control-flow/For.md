# For

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
For[start, test, incr, body] executes start, then repeatedly evaluates body and incr until test fails to give True.
For[start, test, incr] runs with an empty body, useful when incr or test carries the side-effect.
For has attribute HoldAll: start, test, incr, and body are all held unevaluated until For drives them.
Break[] inside body exits the loop.
Continue[] inside body skips the rest of body and re-evaluates incr then test.
Return[v] inside body causes the enclosing function to yield v; For itself returns Null.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Evaluates its arguments in a nonstandard way (sequence: `test`, `body`, `incr`).
- Has attribute `HoldAll`.
- `Break[]` exits the loop.
- `Continue[]` skips the rest of the body and proceeds to evaluating `incr`.
- Exits as soon as `test` fails.
- Returns `Null` unless an explicit `Return` is evaluated.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
