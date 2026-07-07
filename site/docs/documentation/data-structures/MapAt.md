# MapAt

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
MapAt[f, expr, n]
    applies f to the element at position n in expr. Negative n counts from the end.
MapAt[f, expr, {i, j, ...}]
    applies f to the part of expr at position {i, j, ...}.
MapAt[f, expr, {{i1, j1, ...}, {i2, j2, ...}, ...}]
    applies f to the parts of expr at each of the listed positions.

Positions may contain All or Span specifications. MapAt[f, expr, 0]
applies f to the head of expr. Repeated positions apply f repeatedly.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= MapAt[#^2 &, <|"a" -> 3, "b" -> 4|>, "b"]
Out[1]= <|"a" -> 3, "b" -> 16|>

In[2]:= p = <|"a" -> 1, "b" -> 9|>; MapAt[-# &, p, First[Position[p, 9]]]
Out[2]= <|"a" -> 1, "b" -> -9|>
```

## Implementation notes

**Algorithm.** `builtin_map_at` applies `f` at explicit positions rather than at
levels. It first disambiguates a single position from a list of positions: the
argument is treated as *multiple* paths only when it is a non-empty `List` whose
first element is itself a `List`. A single path is then a position vector
`{i1, i2, ...}` (or a bare index), and the recursive `mapat_at_path` walks it:
when the path is exhausted it applies `f` to the targeted node
(`mapat_apply_f`, which builds `f[node]` and calls `evaluate()`); otherwise it
rebuilds the current `EXPR_FUNCTION` with the chosen child replaced by the
recursive result. A path step may be a positive/negative integer (negatives
count from the end, `0` targets the head), the symbol `All` (apply to every
child at that level), or a `Span[a, b]` / `Span[a, b, step]` range. Out-of-range
indices are silently ignored, matching Mathematica's permissive behaviour.

For the multiple-positions form the paths are applied **sequentially** to a
running copy of the expression, so repeated positions apply `f` more than once.

**Data structures.** Operates on the `Expr` tree; copies each level's argument
array and overwrites only the targeted slot, then rebuilds with
`expr_new_function`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= MapAt[f, {a, b, c, d}, 2]
Out[1]= {a, f[b], c, d}
```

```mathematica
In[1]:= MapAt[f, {a, b, c, d}, -1]
Out[1]= {a, b, c, f[d]}
```

```mathematica
In[1]:= MapAt[f, {{a, b}, {c, d}}, {2, 1}]
Out[1]= {{a, b}, {f[c], d}}
```

```mathematica
In[1]:= MapAt[f, {a, b, c, d}, {{1}, {3}}]
Out[1]= {f[a], b, f[c], d}
```

```mathematica
In[1]:= MapAt[Framed, {1, 2, 3, 4, 5}, {{1}, {-1}}]
Out[1]= {Framed[1], 2, 3, 4, Framed[5]}
```

### Notes

`MapAt[f, expr, n]` applies `f` to the element at position `n`; negative `n` counts from the end and `0` targets the head. A position list like `{2, 1}` selects a part deep in a nested expression, while a list of positions `{{1}, {3}}` applies `f` at several places at once (last two examples). Positions may also use `All` or `Span`. Repeated positions apply `f` more than once at that part.
