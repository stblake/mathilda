# Take

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Take[list, n]
    gives the first n elements of list.
Take[list, -n]
    gives the last n elements.
Take[list, {m, n}]
    gives elements m through n.
Take[list, {m, n, s}]
    gives elements m through n in steps of s.
Take[list, {m}]
    gives the single element at position m (wrapped in the head of list).
Take[list, spec1, spec2, ...]
    takes elements at successive levels, e.g. a sub-block of a matrix.

Negative indices count from the end; UpTo[n], All, and None are also accepted
as specifications. Indices are 1-based; out-of-range requests leave the
expression unevaluated. Take operates on any expression, not just List.
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
In[1]:= Take[{a, b, c, d, e}, 3]
Out[1]= {a, b, c}
```

```mathematica
In[1]:= Take[{a, b, c, d, e}, -2]
Out[1]= {d, e}
```

```mathematica
In[1]:= Take[Range[10], {2, 8, 2}]
Out[1]= {2, 4, 6, 8}
```

```mathematica
In[1]:= Take[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 2, 2]
Out[1]= {{1, 2}, {4, 5}}
```

```mathematica
In[1]:= Take[Table[Fibonacci[n], {n, 1, 15}], {3, 15, 3}]
Out[1]= {2, 8, 34, 144, 610}
```

### Notes

`Take[list, n]` takes the first `n` elements, `Take[list, -n]` the last `n`, and
`Take[list, {m, n}]` (optionally `{m, n, s}` with a step) an inclusive index
range. Indices are 1-based and negative indices count from the end; `UpTo[n]`,
`All`, and `None` are also accepted. Multiple specifications act level by level,
so `Take[mat, 2, 2]` extracts the top-left 2x2 sub-block of a matrix. `Take`
operates on any expression, not only `List`; out-of-range requests are left
unevaluated.
