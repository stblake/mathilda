# SameQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
lhs === rhs or SameQ[lhs, rhs]
    yields True if lhs and rhs are structurally identical (head-by-head,
    argument-by-argument), and False otherwise.  Numerically equal but
    distinct heads (e.g. 1 and 1.) are NOT considered same.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_sameq` tests purely structural identity. With fewer than two arguments it returns `True` by convention; otherwise it compares every argument against the first with `expr_eq` and returns `True` only if all are structurally equal, `False` otherwise. Unlike `Equal`, there is no numeric coercion — `1 === 1.` is `False` — and it never returns `NULL` (the result is always a definite boolean).

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/comparisons.c`](https://github.com/stblake/mathilda/blob/main/src/comparisons.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
