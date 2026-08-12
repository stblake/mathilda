# SameQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

lhs === rhs or SameQ\[lhs, rhs\] yields True if lhs and rhs are structurally identical (head-by-head, argument-by-argument), and False otherwise.  Numerically equal but distinct heads (e.g. 1 and 1.) are NOT considered same. SameQ tests structure, not value, so unlike Equal it holds for Indeterminate === Indeterminate.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= a === a
Out[1]= True

In[2]:= 1 === 1.
Out[2]= False

In[3]:= f[x] === f[x]
Out[3]= True
```

## Implementation notes

`builtin_sameq` tests purely structural identity. With fewer than two arguments it returns `True` by convention; otherwise it compares every argument against the first with `expr_eq` and returns `True` only if all are structurally equal, `False` otherwise. Unlike `Equal`, there is no numeric coercion — `1 === 1.` is `False` — and it never returns `NULL` (the result is always a definite boolean).

- Unlike `Equal`, `SameQ` never stays symbolic and does not coerce numeric
  types, so `2 === 2.0` is `False`.
- Also unlike `Equal`, `SameQ` compares structure rather than value, so
  `Indeterminate === Indeterminate` is `True`.

**Attributes:** `Protected`.

## See also

[Equal](../../comparisons/Equal/)

## References

- Source: [`src/comparisons.c`](https://github.com/stblake/mathilda/blob/main/src/comparisons.c)
- Specification: [`docs/spec/builtins/comparisons.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/comparisons.md)
- Tests: [`tests/test_comparisons.c`](https://github.com/stblake/mathilda/blob/main/tests/test_comparisons.c)
- Tests: [`tests/test_hermitian_matrix_q.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hermitian_matrix_q.c)
- Tests: [`tests/test_symmetric_matrix_q.c`](https://github.com/stblake/mathilda/blob/main/tests/test_symmetric_matrix_q.c)

## Notes & additional examples

### Notes

`SameQ` (`===`) is a structural test that always returns `True` or `False`. Numerically equal but distinct heads, such as `1` (Integer) and `1.` (Real), are not the same.
