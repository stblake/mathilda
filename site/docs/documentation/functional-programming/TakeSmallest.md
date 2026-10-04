# TakeSmallest

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TakeSmallest[list, n]`**

Gives the n smallest elements of list, in ascending order. Over an association, gives the n entries with the smallest values (as an association).

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

The three smallest, ascending

```mathematica
In[4]:= TakeSmallest[{3, 1, 4, 1, 5, 9, 2, 6}, 3]
Out[4]= {1, 1, 2}
```

Over an association, the entries with the smallest values

```mathematica
In[5]:= TakeSmallest[<|a -> 3, b -> 1, c -> 5|>, 2]
Out[5]= <|b -> 1, a -> 3|>
```

## Implementation notes

**Algorithm.** `builtin_take_smallest[list, n]` returns the `n` smallest elements
in ascending order — the mirror of `TakeLargest`. For a plain list it calls the
shared `take_extreme(coll, f = NULL, n, largest = false)`: `(key, payload)` pairs
keyed by the element, `qsort`ed ascending (canonical order), then `k = min(n,
len)` payloads copied from the bottom. Over an association the subjects are the
values and the result is the association of the matching entries.

A packed `NDArray` first argument takes the `ndstruct_take_extreme` fast path
(`largest = false`), an `O(n log k)` bounded-heap selection over the raw buffer,
degrading to delist-and-re-evaluate only if that declines.

**Data structures.** `SortByPair` array of `(key, payload)` `Expr` pairs for the
list path; a size-`k` heap over machine numbers for the NDArray path. The input
argument array is borrowed.

**Complexity / limits.** List path `O(n log n)`; NDArray path `O(n log k)`. `n`
must be an explicit integer; a larger `n` returns every element.

**Attributes:** `Protected`.

## References

**See also:** [TakeLargest](../../functional-programming/TakeLargest/), [TakeLargestBy](../../functional-programming/TakeLargestBy/), [TakeSmallestBy](../../functional-programming/TakeSmallestBy/)

- Source: [`src/sort.c`](https://github.com/stblake/mathilda/blob/main/src/sort.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)

## Notes & additional examples

### Notes

`TakeSmallest[list, n]` gives the `n` smallest elements of `list` in ascending
order, the mirror of `TakeLargest`. Over an association it returns the `n` entries
with the smallest values. Use `TakeSmallestBy` to rank by a key function. If `n`
exceeds the length, all elements are returned.
