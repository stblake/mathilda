# ReverseSort

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ReverseSort[list]`**

Sorts into descending order (Reverse of Sort). Over an association, sorts the entries by value, descending.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= ReverseSort[{3, 1, 4, 1, 5, 9, 2}]
Out[1]= {9, 5, 4, 3, 2, 1, 1}

In[2]:= ReverseSort[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[2]= <|"a" -> 3, "c" -> 2, "b" -> 1|>
```

### Applications (3)

Descending order

```mathematica
In[3]:= ReverseSort[{3, 1, 4, 1, 5, 9, 2, 6}]
Out[3]= {9, 6, 5, 4, 3, 2, 1, 1}
```

Canonical order works on strings too

```mathematica
In[4]:= ReverseSort[{"banana", "apple", "cherry"}]
Out[4]= {"cherry", "banana", "apple"}
```

Descending by value, ties kept in input order

```mathematica
In[5]:= ReverseSort[<|a -> 2, b -> 2, c -> 1|>]
Out[5]= <|a -> 2, b -> 2, c -> 1|>
```

## Implementation notes

**Algorithm.** `builtin_reverse_sort` is `Reverse` of `Sort`: it re-heads its
arguments onto `Sort` (`sort_call_as`, so the ascending routine runs as `Sort`
even if it re-evaluates itself), calls `builtin_sort`, and flips the top level with
`reverse_top_level`. The one exception is `ReverseSort[assoc]`, which sorts the
entries descending by value with **equal values kept in input order**
(`sort_collection_p`) — not `Reverse[Sort[assoc]]`, which would flip the ties,
matching Mathematica 15.

`reverse_top_level` handles both list representations. Its `NDArray` arm is a
correctness fix rather than an optimisation: an `NDArray` is `EXPR_NDARRAY`, not
`EXPR_FUNCTION`, so a plain early return left `ReverseSort[NDArray[...]]` showing
`Sort`'s ascending answer; it instead reverses whole rows in place via a
row-sized `memcpy` swap.

**Data structures.** A freshly built `Sort[...]` call (adopting copies of the
arguments), then an in-place pointer reversal of the sorted list's argument array,
or a row-block `memcpy` swap for a packed buffer.

**Complexity / limits.** `O(n log n)` from the underlying `Sort` plus an `O(n)`
reversal. Ordering is Mathilda's canonical order.

**Attributes:** `Protected`.

## References

**See also:** [ReverseSortBy](../../functional-programming/ReverseSortBy/), [Sort](../../data-structures/Sort/)

- Source: [`src/sort.c`](https://github.com/stblake/mathilda/blob/main/src/sort.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_linalg.c)

## Notes & additional examples

### Notes

`ReverseSort[list]` sorts into descending order — `Reverse[Sort[list]]` — using
Mathilda's canonical order, so it ranks numbers, strings, and general expressions
alike. Over an association it sorts the entries by value, descending, keeping
equal values in their input order (so `a -> 2, b -> 2` are not swapped); this is
the Mathematica 15 convention, which is subtly *not* `Reverse[Sort[...]]`.
