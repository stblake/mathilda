# RotateRight

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
RotateRight[expr, n] rotates the elements of expr n positions to the right.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_rotateright` cyclically shifts elements toward the back by `n` (default
1). It negates the shift spec (an integer or a per-level `List` of integers) and delegates to
the same `rotate_rec` worker used by `RotateLeft`, so the offset wrapping and per-level nested
behaviour are identical with the opposite sign.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= RotateRight[{1,2,3,4},1]
Out[1]= {4, 1, 2, 3}
```

### Notes

Cyclically shifts elements `n` positions to the right; it is the inverse of `RotateLeft`.
