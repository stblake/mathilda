# TakeWhile

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TakeWhile[list, crit]`**

Gives the longest leading run of elements e for which crit\[e\] is True. Over an association, tests the values and keeps the matching leading entries (keys preserved).

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= TakeWhile[<|"a" -> 1, "b" -> 2, "c" -> 5, "d" -> 1|>, # < 3 &]
Out[1]= <|"a" -> 1, "b" -> 2|>

In[2]:= LengthWhile[<|"a" -> 1, "b" -> 2, "c" -> 5|>, # < 3 &]
Out[2]= 2
```

### Applications (2)

Stops at 1; the trailing 8 is not taken

```mathematica
In[3]:= TakeWhile[{2, 4, 6, 1, 8}, EvenQ]
Out[3]= {2, 4, 6}
```

```mathematica
In[4]:= TakeWhile[{1, 2, 3, 4}, # < 3 &]
Out[4]= {1, 2}
```

## Implementation notes

**Algorithm.** `builtin_takewhile` returns the leading run of elements satisfying
`crit`, keeping the collection's head. Like `LengthWhile` it first tries the
compiled-predicate fast path (`pred_run_length` finds the run length on the buffer,
then `pred_leading_slice` copies it, so neither the scan nor the result materialises
boxed nodes). A visible `NDArray` is materialised and repacked
(`ndstruct_delist_repack`). Otherwise `leading_run_length` locates the run and the
first `k` elements are copied into a fresh collection; for an association the values
are tested and the head is `Association`.

**Data structures.** A packed buffer on the compiled path; an `Expr**` copy of the
leading `k` elements on the interpreted path.

**Complexity / limits.** O(k) for a run of length `k` (the scan stops at the first
failure), plus O(k) to copy the slice. `LengthWhile` returns only the count.

**Attributes:** `Protected`.

## References

**See also:** [LengthWhile](../../data-structures/LengthWhile/)

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_ndarray_selection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_selection.c)

## Notes & additional examples

### Notes

`TakeWhile` returns the leading run of passing elements and stops at the first
failure — later passing elements are not included. Over an association it tests the
values and returns an association. `LengthWhile` gives the length of the same run.
