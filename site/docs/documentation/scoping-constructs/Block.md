# Block

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Block[{x, y, ...}, expr] evaluates expr with local values for x, y, ....
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `HoldAll`, `Protected`.
- Affects only values, not names.
- Restores original values and attributes after execution.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/scoping-constructs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/scoping-constructs.md)
