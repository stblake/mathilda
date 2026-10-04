# ReverseSortBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ReverseSortBy[list, f]`**

Sorts by f in descending order. Over an association, sorts by f of each value, descending.

**`ReverseSortBy[list, f, p]`**

Sorts by f using the reversed ordering function p.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= ReverseSort[{3, 1, 4, 1, 5, 9, 2}]
Out[1]= {9, 5, 4, 3, 2, 1, 1}

In[2]:= ReverseSort[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[2]= <|"a" -> 3, "c" -> 2, "b" -> 1|>
```

### Applications (2)

Descending by absolute value

```mathematica
In[3]:= ReverseSortBy[{1, -5, 3, -2}, Abs]
Out[3]= {-5, 3, -2, 1}
```

Order the pairs by their last entry, descending

```mathematica
In[4]:= ReverseSortBy[{{1, 2}, {3, 1}, {2, 5}}, Last]
Out[4]= {{2, 5}, {1, 2}, {3, 1}}
```

## Implementation notes

**Algorithm.** `builtin_reverse_sort_by` is `Reverse` of `SortBy`: it re-heads its
arguments onto `SortBy` (`sort_call_as`), calls `builtin_sort_by`, and flips the
top level with the shared `reverse_top_level` helper. The three-argument form
`ReverseSortBy[coll, f, p]` is treated as `SortBy[coll, f, p[#2, #1] &]`
(`sort_by_with_p` with the ordering reversed), which differs from
`Reverse[SortBy[coll, f, p]]` on ties — again the Mathematica-15 convention.

As with `ReverseSort`, the `NDArray` arm of `reverse_top_level` is a correctness
fix: a packed buffer is reversed by a row-sized `memcpy` swap rather than silently
passing through as `SortBy`'s ascending result. Over an association the sort is by
`f` of each value.

**Data structures.** A synthesised `SortBy[...]` call (copies of the arguments)
and an in-place reversal of the sorted result's argument array (or a row-block
`memcpy` swap for a packed buffer).

**Complexity / limits.** `O(n)` key evaluations plus the `O(n log n)` underlying
`SortBy` and an `O(n)` reversal.

**Attributes:** `Protected`.

## References

**See also:** [ReverseSort](../../functional-programming/ReverseSort/), [Sort](../../data-structures/Sort/)

- Source: [`src/sort.c`](https://github.com/stblake/mathilda/blob/main/src/sort.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

`ReverseSortBy[list, f]` sorts by `f` in descending order — the `By` companion to
`ReverseSort`, and the descending counterpart of `SortBy`. The key function `f` is
applied to each element (or each value, for an association) and the elements are
ordered by the canonical order of those keys, largest first. A three-argument
`ReverseSortBy[list, f, p]` ranks the keys with the reversed ordering function `p`.
