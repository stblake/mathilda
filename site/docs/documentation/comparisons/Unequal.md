# Unequal

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

lhs != rhs or Unequal\[lhs, rhs\] is the negation of Equal: True if lhs and rhs can be decided unequal, False if they can be decided equal, otherwise unevaluated. Per IEEE 754, an Indeterminate argument gives True.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= 2 != 3
Out[1]= True

In[2]:= 2 != 2
Out[2]= False

In[3]:= a != b
Out[3]= a != b
```

## Implementation notes

`builtin_unequal` is the value-level (not structural) negation of `Equal`. For every argument pair it first checks `expr_eq`; on failure it tries `compare_numeric` to decide equality/inequality numerically. If any pair is found equal it returns `False` immediately. It returns `True` only when *every* pair is provably unequal — either decided by `compare_numeric` or, for non-comparable values, when both sides are distinct raw data (`is_raw_data`). If some pair is neither equal nor provably unequal (e.g. symbolic), it returns `NULL` so the call stays unevaluated. Fewer than two arguments returns `True`.

- A pair containing `Indeterminate` counts as unequal, so
  `Indeterminate != Indeterminate` is `True`. An equal non-`Indeterminate` pair
  still decides the whole call `False`.
- Like `Equal`, an undecided pair of numeric constants is settled by the exact
  zero-test, so `I != 0` and `2 I Pi != 0` are `True`, while a free symbol stays
  symbolic (`x != 0`).

**Attributes:** `Protected`.

## References

**See also:** [Equal](../../comparisons/Equal/)

- Source: [`src/comparisons.c`](https://github.com/stblake/mathilda/blob/main/src/comparisons.c)
- Specification: [`docs/spec/builtins/comparisons.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/comparisons.md)
- Tests: [`tests/test_comparisons.c`](https://github.com/stblake/mathilda/blob/main/tests/test_comparisons.c)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)

## Notes & additional examples

### Notes

`Unequal` (`!=`) is the negation of `Equal`: `True` if the arguments are decidably unequal, `False` if decidably equal, otherwise it stays unevaluated.
