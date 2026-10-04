# TakeLargestBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TakeLargestBy[list, f, n]`**

Gives the n elements of list for which f is largest, in descending order of f. Over an association, ranks by f of each value.

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

The two with the largest absolute value

```mathematica
In[4]:= TakeLargestBy[{1, -5, 3, -2, 4}, Abs, 2]
Out[4]= {-5, 4}
```

Rank the pairs by their last entry

```mathematica
In[5]:= TakeLargestBy[{{1, 9}, {5, 2}, {3, 7}}, Last, 2]
Out[5]= {{1, 9}, {3, 7}}
```

## Implementation notes

**Algorithm.** `builtin_take_largest_by[list, f, n]` returns the `n` elements for
which `f` is largest, in descending order of `f`. It calls the shared
`take_extreme(coll, f, n, largest = true)`: a `(key, payload)` pair is formed per
element with the key `f[subject]` evaluated eagerly (for an association, `f` of
each value), the pairs are `qsort`ed ascending by key in canonical order, and the
top `k = min(n, len)` payloads are copied walking down from the largest key. The
original elements — not the keys — are returned, keeping the collection's head.

**Data structures.** `SortByPair` array of `n` `(key, payload)` `Expr` pairs; each
key is a freshly built and evaluated `f[subject]`. All pairs are freed once the
result node is built; the input argument array is borrowed.

**Complexity / limits.** `O(n)` key evaluations plus an `O(n log n)` sort (the
`By` form does not take the bounded-heap NDArray fast path that keyless
`TakeLargest` uses). `n` must be an explicit integer; a larger `n` returns every
element.

**Attributes:** `Protected`.

## References

**See also:** [TakeLargest](../../functional-programming/TakeLargest/), [TakeSmallest](../../functional-programming/TakeSmallest/), [TakeSmallestBy](../../functional-programming/TakeSmallestBy/)

- Source: [`src/sort.c`](https://github.com/stblake/mathilda/blob/main/src/sort.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

`TakeLargestBy[list, f, n]` gives the `n` elements of `list` for which `f` is
largest, in descending order of `f`. The ranking uses `f` of each element but the
original elements are returned — here `-5` and `4` have the two largest absolute
values. Over an association it ranks by `f` of each value. If `n` exceeds the
length, all elements are returned.
