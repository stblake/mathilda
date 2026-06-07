# Array

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Array[f, n]
    generates a list {f[1], f[2], ..., f[n]}.
Array[f, n, r]
    generates a list of length n starting from index r.
Array[f, {n1, n2, ...}]
    generates an n1 x n2 x ... nested-list array with elements
    f[i1, i2, ...].
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
