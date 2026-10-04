# DeleteDuplicatesBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DeleteDuplicatesBy[expr, f]`**

Keeps the first element for each distinct f\[element\], preserving order. Over an association, f is applied to the values and the surviving entries are kept (keys preserved).

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= DeleteDuplicatesBy[{1, 12, 3, 14, 5}, EvenQ]
Out[1]= {1, 12}

In[2]:= DeleteDuplicatesBy[<|"a" -> 1, "b" -> 12, "c" -> 3|>, EvenQ]
Out[2]= <|"a" -> 1, "b" -> 12|>
```

### Applications (2)

First odd and first even survive

```mathematica
In[3]:= DeleteDuplicatesBy[{1, 2, 3, 4, 5, 6}, EvenQ]
Out[3]= {1, 2}
```

Keyed by |x|

```mathematica
In[4]:= DeleteDuplicatesBy[{-1, 1, 2, -2, 3}, Abs]
Out[4]= {-1, 2, 3}
```

## Implementation notes

**Algorithm.** `builtin_deleteduplicatesby` keeps the first element for each
distinct value of `f[element]`, preserving order. It walks the collection,
evaluates `f` on each element (on each *value*, for an association), and compares
the result against the `f`-values of the survivors so far; a new `f`-value adds
the element to the kept list, a repeat drops it. Over an association the surviving
entries are returned as an association with keys preserved.

**Data structures.** Two parallel arrays — the surviving elements/rules and their
owned `f`-values — with survivors compared directly by `expr_eq`. The survivor
count is typically small, so the linear scan of seen `f`-values stays within
budget; no hash index is built.

**Complexity / limits.** `O(n · k)` for `n` elements and `k` survivors
(`O(n²)` worst case when nearly everything is distinct), plus `n` evaluations of
`f`. Requires a 2-argument call with a non-atomic first argument; `f` is
arbitrary, so no packed fast path applies.

**Attributes:** `Protected`.

## References

- Source: [`src/list/setops.c`](https://github.com/stblake/mathilda/blob/main/src/list/setops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

`DeleteDuplicatesBy[list, f]` keeps the first element for each distinct value of
`f[element]`, preserving order — the deduplicating cousin of `GatherBy`, returning
one representative per group instead of the groups themselves. With
`f = Abs` it collapses `±x` to whichever sign appears first. Over an association
`f` is applied to each value and the surviving entries are returned as an
association with keys preserved. Use `DeleteDuplicates` for the plain
`f = Identity` case.
