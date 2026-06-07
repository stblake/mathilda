# ListQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ListQ[expr] gives True if expr is a list (head List), False otherwise.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_listq` (`src/list.c`) returns `True`/`False` according to the `is_listq` helper, i.e. whether the argument is an `EXPR_FUNCTION` whose head is the symbol `List`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
