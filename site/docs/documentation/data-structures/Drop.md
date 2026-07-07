# Drop

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Drop[list, n]
    gives list with its first n elements dropped.
Drop[list, -n]
    drops the last n elements.
Drop[list, {m, n}]
    drops elements m through n.
Drop[list, {m, n, s}]
    drops elements m through n in steps of s.
Drop[list, {m}]
    drops the single element at position m.
Drop[list, spec1, spec2, ...]
    drops elements at successive levels.

Negative indices count from the end; UpTo[n], All, and None are also accepted.
Indices are 1-based; out-of-range requests leave the expression unevaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= First[<|"a" -> 10, "b" -> 20|>]
Out[1]= 10

In[2]:= Rest[<|"a" -> 10, "b" -> 20, "c" -> 30|>]
Out[2]= <|"b" -> 20, "c" -> 30|>

In[3]:= Take[<|"a" -> 1, "b" -> 2, "c" -> 3|>, 2]
Out[3]= <|"a" -> 1, "b" -> 2|>
```

## Implementation notes

**Attributes:** `NHoldRest`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)

## Notes & additional examples

### Worked examples

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

### Notes

`Drop` is the complement of `Take`. A plain count drops from the front
(`Drop[list, n]`) or, with a negative count, from the back. The `{m, n}` form
drops a contiguous block, `{m, n, s}` drops a strided slice, and `{m}` drops a
single element. Multiple level specifications drop along successive list
dimensions, so the `{2}, {2}` example deletes the second row and the second
column of a matrix in one call. Indices are 1-based and negative indices count
from the end; out-of-range requests leave the expression unevaluated.
