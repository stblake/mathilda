# UnsameQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
lhs =!= rhs or UnsameQ[lhs, rhs]
    is the negation of SameQ: True iff lhs and rhs are not structurally
    identical.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_unsameq` is the negation of structural identity over all argument pairs. With fewer than two arguments it returns `True`; otherwise it checks every pair `(i,j)` with `expr_eq` and returns `False` as soon as any two are structurally equal, else `True`. Like `SameQ` it uses no numeric coercion and always yields a definite boolean.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/comparisons.c`](https://github.com/stblake/mathilda/blob/main/src/comparisons.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
