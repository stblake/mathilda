# TakeLargest

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TakeLargest[list, n]`**

Gives the n largest elements of list, in descending order. Over an association, gives the n entries with the largest values (as an association).

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

The three largest, descending

```mathematica
In[4]:= TakeLargest[{3, 1, 4, 1, 5, 9, 2, 6}, 3]
Out[4]= {9, 6, 5}
```

Over an association, the entries with the largest values

```mathematica
In[5]:= TakeLargest[<|a -> 3, b -> 1, c -> 5|>, 2]
Out[5]= <|c -> 5, a -> 3|>
```

## Implementation notes

**Algorithm.** `builtin_take_largest[list, n]` returns the `n` largest elements in
descending order. For a plain list it calls the shared `take_extreme(coll, f =
NULL, n, largest = true)`: it forms `(key, payload)` pairs with the key being the
element itself, `qsort`s them ascending by key (`sortby_pair_cmp`, canonical
order), and copies `k = min(n, len)` payloads walking down from the top. Over an
association the subjects are the values and the result is the association of the
matching entries.

When the first argument is a packed `NDArray`, it takes the `ndstruct_take_extreme`
fast path instead — an `O(n log k)` bounded-heap selection over the raw buffer, with
no boxed `Expr` per element — and degrades to delist-and-re-evaluate only if that
declines.

**Data structures.** `SortByPair` array of `n` `(key, payload)` `Expr` pairs for
the list path (freed after the result is built); a size-`k` heap over machine
numbers for the NDArray path. The input argument array is borrowed.

**Complexity / limits.** List path `O(n log n)` (a full sort). NDArray path
`O(n log k)` through the heap. `n` must be an explicit integer; if it exceeds the
length, all elements are returned.

**Attributes:** `Protected`.

## References

**See also:** [TakeSmallest](../../functional-programming/TakeSmallest/), [TakeLargestBy](../../functional-programming/TakeLargestBy/), [TakeSmallestBy](../../functional-programming/TakeSmallestBy/)

- Source: [`src/sort.c`](https://github.com/stblake/mathilda/blob/main/src/sort.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)

## Notes & additional examples

### Notes

`TakeLargest[list, n]` gives the `n` largest elements of `list` in descending
order, by Mathilda's canonical order. Over an association it returns the `n`
entries with the largest values, as an association. Use `TakeLargestBy` to rank by
a key function instead of by the elements themselves. If `n` exceeds the length,
all elements are returned.
