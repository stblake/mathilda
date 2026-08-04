# FlattenAt

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FlattenAt[list, n]
    flattens out the sublist at position n of list, splicing its elements
    into list; a negative n counts from the end.
FlattenAt[expr, {i, j, ...}]
    flattens out the part of expr at the position {i, j, ...}.
FlattenAt[expr, {{i1, ...}, {i2, ...}, ...}]
    flattens out the parts of expr at several positions.
    The head of the spliced part is removed; FlattenAt works on any head,
    not just List.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
