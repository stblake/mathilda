# SplitBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SplitBy[list, f]
    splits list into runs of consecutive elements that give the same
    value of f[element]. Only adjacent elements are grouped (unlike
    GatherBy, which collects equal keys from anywhere in the list).
SplitBy[list, {f1, f2, ...}]
    splits by f1, then splits each resulting run by f2, and so on,
    nesting one level deeper per function.
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
