# RotateRight

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RotateRight[expr, n] rotates the elements of expr n positions to the right.`**

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (1)

```mathematica
In[1]:= RotateRight[{1,2,3,4},1]
Out[1]= {4, 1, 2, 3}
```

## Implementation notes

**Algorithm.** `builtin_rotateright` cyclically shifts elements toward the back by `n` (default
1). It negates the shift spec (an integer or a per-level `List` of integers) and delegates to
the same `rotate_rec` worker used by `RotateLeft`, so the offset wrapping and per-level nested
behaviour are identical with the opposite sign.

**Attributes:** `Protected`.

## See also

[RotateLeft](../../structural-manipulation/RotateLeft/)

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_ndarray_selection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_selection.c)

## Notes & additional examples

### Notes

Cyclically shifts elements `n` positions to the right; it is the inverse of `RotateLeft`.
