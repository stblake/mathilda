# Drop

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Drop[list, n]`**

gives list with its first n elements dropped.

**`Drop[list, -n]`**

drops the last n elements.

**`Drop[list, {m, n}]`**

drops elements m through n.

**`Drop[list, {m, n, s}]`**

drops elements m through n in steps of s.

**`Drop[list, {m}]`**

drops the single element at position m.

**`Drop[list, spec1, spec2, ...]`**

drops elements at successive levels.

<details>
<summary>Notes</summary>

Negative indices count from the end; UpTo\[n\], All, and None are also accepted. Indices are 1-based; out-of-range requests leave the expression unevaluated.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Rest[<|"a" -> 10, "b" -> 20, "c" -> 30|>]
Out[1]= <|"b" -> 20, "c" -> 30|>

In[2]:= Take[<|"a" -> 1, "b" -> 2, "c" -> 3|>, 2]
Out[2]= <|"a" -> 1, "b" -> 2|>
```

### Applications (5)

```mathematica
In[1]:= Drop[{a, b, c, d, e}, 2]
Out[1]= {c, d, e}
```

```mathematica
In[1]:= Drop[{a, b, c, d, e}, -2]
Out[1]= {a, b, c}
```

```mathematica
In[1]:= Drop[{a, b, c, d, e}, {2, 4}]
Out[1]= {a, e}
```

```mathematica
In[1]:= Drop[{a, b, c, d, e, f, g}, {2, 7, 2}]
Out[1]= {a, c, e, g}
```

```mathematica
In[1]:= Drop[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {2}, {2}]
Out[1]= {{1, 3}, {7, 9}}
```

## Implementation notes

**Attributes:** `NHoldRest`, `Protected`.

## See also

[First](../../data-structures/First/), [Last](../../data-structures/Last/), [Rest](../../data-structures/Rest/), [Most](../../data-structures/Most/), [Take](../../data-structures/Take/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_ndarray_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_reduce.c)

## Notes & additional examples

### Notes

`Drop` is the complement of `Take`. A plain count drops from the front
(`Drop[list, n]`) or, with a negative count, from the back. The `{m, n}` form
drops a contiguous block, `{m, n, s}` drops a strided slice, and `{m}` drops a
single element. Multiple level specifications drop along successive list
dimensions, so the `{2}, {2}` example deletes the second row and the second
column of a matrix in one call. Indices are 1-based and negative indices count
from the end; out-of-range requests leave the expression unevaluated.
