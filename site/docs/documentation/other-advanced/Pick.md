# Pick

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Pick[expr, sel]
    Picks out the elements of expr for which the
    corresponding element of sel is True.
Pick[expr, sel, patt]
    Picks out the elements of expr for which the
    corresponding element of sel matches patt.
    Operates at all levels; sel must mirror the structure of expr, and
    the head of expr is preserved. Returns unevaluated if the structures
    disagree.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
