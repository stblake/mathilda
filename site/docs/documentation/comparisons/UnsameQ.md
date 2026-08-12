# UnsameQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

lhs =!= rhs or UnsameQ\[lhs, rhs\] is the negation of SameQ: True iff lhs and rhs are not structurally identical.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

```mathematica
In[1]:= a =!= b
Out[1]= True

In[2]:= a =!= a
Out[2]= False
```

## Implementation notes

`builtin_unsameq` is the negation of structural identity over all argument pairs. With fewer than two arguments it returns `True`; otherwise it checks every pair `(i,j)` with `expr_eq` and returns `False` as soon as any two are structurally equal, else `True`. Like `SameQ` it uses no numeric coercion and always yields a definite boolean.

**Attributes:** `Protected`.

## See also

[SameQ](../../comparisons/SameQ/)

## References

- Source: [`src/comparisons.c`](https://github.com/stblake/mathilda/blob/main/src/comparisons.c)
- Specification: [`docs/spec/builtins/comparisons.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/comparisons.md)
- Tests: [`tests/test_comparisons.c`](https://github.com/stblake/mathilda/blob/main/tests/test_comparisons.c)

## Notes & additional examples

### Notes

`UnsameQ` (`=!=`) is the structural negation of `SameQ`; it always decides to `True` or `False`.
