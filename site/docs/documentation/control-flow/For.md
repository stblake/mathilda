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

**Algorithm.** `builtin_for` implements the C-style `For[start, test, incr]` or `For[start, test, incr, body]` and is `ATTR_HOLDALL`, so each part is re-evaluated every pass. It evaluates `start` once for its side effects (discarding the result), then loops: `evaluate(test)`; exit unless it is exactly the symbol `True`; `evaluate(body)`; `evaluate(incr)`. Both the body and the increment are passed through `iter_flow_classify` with boundary head `SYM_For`: a `Break` exits (loop yields `Null`), a `Return`/targeted-Return payload exits with that value, and `Throw`/`Abort`/`Quit` propagate unchanged. `Continue` inside the body falls through to the increment step; `Continue` inside the increment simply proceeds to the next test. Iterator variables are not localized (the user mutates ordinary symbols in `start`/`incr`). Returns the Return payload if one was issued, else `Null`.

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

- Source: [`src/iter.c`](https://github.com/stblake/mathilda/blob/main/src/iter.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
