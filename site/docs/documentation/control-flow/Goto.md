# Goto

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Goto[tag]
    Scans the CompoundExpression it appears in directly for
    Label[tag], then enclosing ones, and transfers control to that point.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Both are `Protected`. `tag` is evaluated (conventionally a literal symbol or

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
