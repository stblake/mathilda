# LengthWhile

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LengthWhile[list, crit]`**

Gives the length of the longest leading run of elements e for which crit\[e\] is True. Over an association, tests values.

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

Leading run stops at -1

```mathematica
In[3]:= LengthWhile[{1, 2, 3, -1, 5}, Positive]
Out[3]= 3
```

Stops at 7, even though 8 follows

```mathematica
In[4]:= LengthWhile[{2, 4, 6, 7, 8}, EvenQ]
Out[4]= 3
```

## Implementation notes

**Algorithm.** `builtin_lengthwhile` returns the length of the leading run of
elements for which `crit` gives `True`, stopping at the first failure. It tries a
compiled-predicate fast path first: `pred_run_length` compiles `crit` and scans a
machine buffer directly, so for a numeric predicate nothing is boxed and the answer
is just the count. Failing that, a visible `NDArray` is materialised and repacked
(`ndstruct_delist_repack`), and otherwise `leading_run_length` walks the elements one
by one (for an association it tests the *values*).

**Data structures.** A packed numeric buffer on the compiled path; plain `Expr`
element pointers on the interpreted path. Only a count is returned, so nothing is
copied.

**Complexity / limits.** O(k) where `k` is the run length — the scan stops at the
first element that fails `crit`, so a predicate that fails immediately is O(1).
`TakeWhile` returns the run itself rather than its length.

**Attributes:** `Protected`.

## References

**See also:** [TakeWhile](../../data-structures/TakeWhile/)

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_ndarray_selection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_selection.c)
- Tests: [`tests/test_pred_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_pred_compile.c)

## Notes & additional examples

### Notes

Only the *leading* run counts: the scan halts at the first element failing the
predicate, so later passing elements are not counted. A numeric predicate takes a
compiled buffer scan that builds nothing. `TakeWhile` returns that leading run as a
collection instead of its length.
