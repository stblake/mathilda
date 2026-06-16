# Flatten

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Flatten[list]
    flattens out nested lists, collapsing every level into a flat list
    with the same head as the top level.
Flatten[list, n]
    flattens only the top n levels.
Flatten[list, n, h]
    flattens only sublists whose head matches h, leaving other heads
    in place.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_flatten` (in `src/list.c`) accepts `Flatten[list]`, `Flatten[list, n]` (level cap), and `Flatten[list, n, h]` (custom head). It iterates over the top-level arguments calling the recursive worker `flatten_rec`, which splices the children of any subexpression whose head equals the flattening head `h` (default `List`) up into the output, descending up to `n` levels (n = -1 means unlimited). The collected arguments are gathered into a growable buffer and reassembled under the original head.

**Data structures.** A dynamically grown `Expr**` accumulator (`results`, with `count`/`cap`) holds the deep-copied leaf expressions before the final `expr_new_function` rebuild.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

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

### Notes

Without a level argument, `Flatten` removes all nesting; with a level `n` it
flattens only the top `n` levels. The third example collapses a 2x2 nested
`Table` into a flat list of coordinate pairs (a common reshaping idiom). The
last shows the three-argument form `Flatten[expr, n, h]`, which flattens
nested calls of an arbitrary head `h` rather than `List`.
