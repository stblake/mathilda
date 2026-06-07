# InputForm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
InputForm[expr]
    prints expr in a form suitable to be re-read by the parser, using
    operator syntax (a + b, not Plus[a, b]) and explicit string quotes.
Like FullForm, InputForm is a printer wrapper: it is consumed during
output and does not appear in the printed result.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`InputForm` is an unevaluated display wrapper: `builtin_inputform` (`src/print.c`) returns `NULL`, leaving `InputForm[expr]` intact. The printer's `print_standard` detects the `InputForm` head and renders the argument in a re-parseable form (a printer flag toggles InputForm-specific formatting); `ToString[expr, InputForm]` routes through the same standard printer.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/print.c`](https://github.com/stblake/mathilda/blob/main/src/print.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
