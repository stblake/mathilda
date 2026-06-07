# Union

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Union[list]
    gives the sorted list of distinct elements in list.
Union[l1, l2, ...]
    gives the sorted list of distinct elements appearing in any of the
    input lists (set union).
Comparison is by canonical structural equality.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Union[{1, 2, 1, 3, 6, 2, 2}]
Out[1]= {1, 2, 3, 6}

In[2]:= Union[{a, b, a, c}, {d, a, e, b}, {c, a}]
Out[2]= {a, b, c, d, e}
```

## Implementation notes

**Algorithm.** `builtin_union` concatenates the elements of all argument lists (which must
share a common head), sorts the combined `Expr**` array with `qsort` under the canonical
`expr_compare` order, then removes adjacent duplicates — `expr_eq` by default, or an optional
`SameTest -> f` which is evaluated per adjacent pair. The result is the sorted, deduplicated
list. (`DeleteDuplicates` in the same file does the order-preserving variant using a hash table
keyed on `expr_hash`/`expr_eq`.)

- `Flat`, `OneIdentity`, `Protected`, `ReadProtected`.
- All expressions must have the same head.
- Result has the same head as the inputs.

**Attributes:** `Flat`, `OneIdentity`, `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
