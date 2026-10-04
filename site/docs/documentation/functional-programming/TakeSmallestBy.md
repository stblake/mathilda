# TakeSmallestBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TakeSmallestBy[list, f, n]`**

Gives the n elements of list for which f is smallest, in ascending order of f. Over an association, ranks by f of each value.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= TakeLargest[{3, 1, 4, 1, 5, 9, 2, 6}, 3]
Out[1]= {9, 6, 5}

In[2]:= TakeLargest[<|"a" -> 3, "b" -> 9, "c" -> 1, "d" -> 6|>, 2]
Out[2]= <|"b" -> 9, "d" -> 6|>

In[3]:= TakeLargestBy[{-9, 2, -3, 5}, Abs, 2]
Out[3]= {-9, 5}
```

### Applications (2)

The two with the smallest absolute value

```mathematica
In[4]:= TakeSmallestBy[{1, -5, 3, -2, 4}, Abs, 2]
Out[4]= {1, -2}
```

Rank the pairs by their last entry

```mathematica
In[5]:= TakeSmallestBy[{{1, 9}, {5, 2}, {3, 7}}, Last, 2]
Out[5]= {{5, 2}, {3, 7}}
```

## Implementation notes

**Algorithm.** `builtin_take_smallest_by[list, f, n]` returns the `n` elements for
which `f` is smallest, in ascending order of `f` — the mirror of `TakeLargestBy`.
It calls the shared `take_extreme(coll, f, n, largest = false)`: a `(key,
payload)` pair per element keyed by the eagerly evaluated `f[subject]` (for an
association, `f` of each value), `qsort`ed ascending in canonical order, and the
bottom `k = min(n, len)` payloads copied from the smallest key up. The original
elements are returned, keeping the collection's head.

**Data structures.** `SortByPair` array of `n` `(key, payload)` `Expr` pairs; each
key is a freshly built and evaluated `f[subject]`. The pairs are freed after the
result is built; the input argument array is borrowed.

**Complexity / limits.** `O(n)` key evaluations plus an `O(n log n)` sort. `n`
must be an explicit integer; a larger `n` returns every element.

**Attributes:** `Protected`.

## References

**See also:** [TakeLargest](../../functional-programming/TakeLargest/), [TakeSmallest](../../functional-programming/TakeSmallest/), [TakeLargestBy](../../functional-programming/TakeLargestBy/)

- Source: [`src/sort.c`](https://github.com/stblake/mathilda/blob/main/src/sort.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)

## Notes & additional examples

### Notes

`TakeSmallestBy[list, f, n]` gives the `n` elements of `list` for which `f` is
smallest, in ascending order of `f` — the mirror of `TakeLargestBy`. The ranking
uses `f` of each element but returns the original elements. Over an association it
ranks by `f` of each value. If `n` exceeds the length, all elements are returned.
