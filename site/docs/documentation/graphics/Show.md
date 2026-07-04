# Show

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Show[graphics, opts...]
    Displays graphics (a Graphics[...] object) in an interactive window and returns it, merging any given options into its option list.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`.
- Declines to evaluate (stays unevaluated) if its first argument isn't a

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
