# Flatten

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Flatten[list]`**

flattens out nested lists, collapsing every level into a flat list with the same head as the top level.

**`Flatten[list, n]`**

flattens only the top n levels.

**`Flatten[list, n, h]`**

flattens only sublists whose head matches h, leaving other heads in place.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= Flatten[{{1, 2}, {3, {4}}}]
Out[1]= {1, 2, 3, 4}
```

```mathematica
In[1]:= Flatten[{{1, {2, 3}}, {4, {5}}}, 1]
Out[1]= {1, {2, 3}, 4, {5}}
```

```mathematica
In[1]:= Flatten[Table[{i, j}, {i, 2}, {j, 2}], 1]
Out[1]= {{1, 1}, {1, 2}, {2, 1}, {2, 2}}
```

```mathematica
In[1]:= Flatten[f[a, f[b, f[c, d]]], 2, f]
Out[1]= f[a, b, c, d]
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Dot 6x6 x 6x6 x 10000 | 338 s | 6.36 s | 4.01 s |
| Eigenvalues 600x600 (general) | 270 s | 163 s | 79.2 s |
| Transpose then Dot (fused?) | 37.1 s | 63.1 s | 18.6 s |
| SingularValueDecomposition 400x400 | 23.6 s | 45.7 s | 8.93 s |
| QRDecomposition 400x400 | 7.98 s | 18.6 s | 2.46 s |
| Inverse 3x3 x 5000 | 2.7 s | 4.75 s | 7.67 s |

## Implementation notes

**Algorithm.** `builtin_flatten` (in `src/list.c`) accepts `Flatten[list]`, `Flatten[list, n]` (level cap), and `Flatten[list, n, h]` (custom head). It iterates over the top-level arguments calling the recursive worker `flatten_rec`, which splices the children of any subexpression whose head equals the flattening head `h` (default `List`) up into the output, descending up to `n` levels (n = -1 means unlimited). The collected arguments are gathered into a growable buffer and reassembled under the original head.

**Data structures.** A dynamically grown `Expr**` accumulator (`results`, with `count`/`cap`) holds the deep-copied leaf expressions before the final `expr_new_function` rebuild.

**Attributes:** `Protected`.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_blas.c`](https://github.com/stblake/mathilda/blob/main/tests/test_blas.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)

## Notes & additional examples

### Notes

Without a level argument, `Flatten` removes all nesting; with a level `n` it
flattens only the top `n` levels. The third example collapses a 2x2 nested
`Table` into a flat list of coordinate pairs (a common reshaping idiom). The
last shows the three-argument form `Flatten[expr, n, h]`, which flattens
nested calls of an arbitrary head `h` rather than `List`.
