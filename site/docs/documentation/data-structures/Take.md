# Take

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Take[list, n]`**

gives the first n elements of list.

**`Take[list, -n]`**

gives the last n elements.

**`Take[list, {m, n}]`**

gives elements m through n.

**`Take[list, {m, n, s}]`**

gives elements m through n in steps of s.

**`Take[list, {m}]`**

gives the single element at position m (wrapped in the head of list).

**`Take[list, spec1, spec2, ...]`**

takes elements at successive levels, e.g. a sub-block of a matrix.

<details>
<summary>Notes</summary>

Negative indices count from the end; UpTo\[n\], All, and None are also accepted as specifications. Indices are 1-based; out-of-range requests leave the expression unevaluated. Take operates on any expression, not just List.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= First[<|"a" -> 10, "b" -> 20|>] In[1b]:= First[<||>, 0] Out[1b]= 0

In[2]:= Rest[<|"a" -> 10, "b" -> 20, "c" -> 30|>]
Out[2]= <|"b" -> 20, "c" -> 30|>

In[3]:= Take[<|"a" -> 1, "b" -> 2, "c" -> 3|>, 2]
Out[3]= <|"a" -> 1, "b" -> 2|>
```

### Applications (5)

```mathematica
In[4]:= Take[{a, b, c, d, e}, 3]
Out[4]= {a, b, c}

In[5]:= Take[{a, b, c, d, e}, -2]
Out[5]= {d, e}

In[6]:= Take[Range[10], {2, 8, 2}]
Out[6]= {2, 4, 6, 8}

In[7]:= Take[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 2, 2]
Out[7]= {{1, 2}, {4, 5}}

In[8]:= Take[Table[Fibonacci[n], {n, 1, 15}], {3, 15, 3}]
Out[8]= {2, 8, 34, 144, 610}
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Sort 4x10^6 | 42.2 s | 68.7 s | 111 s |
| gather v[[idx]], 4x10^6 | 16.8 s | 6.66 s | 7.18 s |
| Union of 4x10^6 integers | 12.4 s | 71.1 s | 376 s |
| Reverse 4x10^6 | 5.37 s | 0.297 s | 0.982 s |
| Join two 2x10^6 | 0.899 s | 0.6 s | 0.397 s |
| RotateLeft 4x10^6 by 1000 | 0.855 s | 0.307 s | 0.456 s |

## Implementation notes

**Attributes:** `NHoldRest`, `Protected`.

## References

**See also:** [First](../../data-structures/First/), [Last](../../data-structures/Last/), [Rest](../../data-structures/Rest/), [Most](../../data-structures/Most/), [Drop](../../data-structures/Drop/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_clip.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clip.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`Take[list, n]` takes the first `n` elements, `Take[list, -n]` the last `n`, and
`Take[list, {m, n}]` (optionally `{m, n, s}` with a step) an inclusive index
range. Indices are 1-based and negative indices count from the end; `UpTo[n]`,
`All`, and `None` are also accepted. Multiple specifications act level by level,
so `Take[mat, 2, 2]` extracts the top-left 2x2 sub-block of a matrix. `Take`
operates on any expression, not only `List`; out-of-range requests are left
unevaluated.
