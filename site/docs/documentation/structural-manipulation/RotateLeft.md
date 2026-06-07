# RotateLeft

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
RotateLeft[expr, n] rotates the elements of expr n positions to the left.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_rotateleft` cyclically shifts elements toward the front by `n`
(default 1) using `rotate_rec`: at each level it computes the wrapped offset `((n mod len) +
len) mod len` and reads element `i` from source index `(i + offset) mod len`. A `List`-valued
`n` applies a per-level shift amount as the recursion descends into nested lists. `RotateRight`
is implemented by negating `n` and calling the same routine.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
