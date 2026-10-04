# FlattenAt

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FlattenAt[list, n]`**

flattens out the sublist at position n of list, splicing its elements into list; a negative n counts from the end.

**`FlattenAt[expr, {i, j, ...}]`**

flattens out the part of expr at the position {i, j, ...}.

**`FlattenAt[expr, {{i1, ...}, {i2, ...}, ...}]`**

flattens out the parts of expr at several positions. The head of the spliced part is removed; FlattenAt works on any head, not just List.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Worked examples (2)

```mathematica
In[1]:= FlattenAt[{a,{b,c},{d,e},{f}}, 2]
Out[1]= {a, b, c, {d, e}, {f}}

In[2]:= FlattenAt[{a,{b,c},{d,e},{f}}, {{2},{4}}]
Out[2]= {a, b, c, {d, e}, f}
```

### Applications (4)

Splice the sublist at position 2 into its parent

```mathematica
In[3]:= FlattenAt[{a, {b, c}, {d, e}, {f}}, 2]
Out[3]= {a, b, c, {d, e}, {f}}
```

Several positions, each against the original

```mathematica
In[4]:= FlattenAt[{a, {b, c}, {d, e}, {f}}, {{2}, {4}}]
Out[4]= {a, b, c, {d, e}, f}
```

Any head works, not just List

```mathematica
In[5]:= FlattenAt[f[g[1, 2], g[3, 4]], 1]
Out[5]= f[1, 2, g[3, 4]]
```

A negative position counts from the end

```mathematica
In[6]:= FlattenAt[{a, {b, c}, {d, e}, {f}}, -1]
Out[6]= {a, {b, c}, {d, e}, f}
```

## Implementation notes

**Algorithm.** `builtin_flatten_at` splices the subexpression at one or more
positions into its parent, removing that subexpression's head — it works on any
head, not just `List` (`FlattenAt[f[g[1, 2], g[3, 4]], 1]` is `f[1, 2, g[3, 4]]`).
Position resolution — an integer (negatives counting from the end), a single deep
path `{i, j, ...}`, or a list of paths `{{i1, ...}, {i2, ...}}` — is delegated to
the shared walker `expr_apply_at_positions` (`src/part.c`), the same one `MapAt`
and `ReplaceAt` use, so the one-path/many-paths distinction is exactly theirs; an
out-of-range position makes the walker return `NULL` and `FlattenAt` stay
unevaluated. The splice itself is deferred to the evaluator: the leaf action
`flatten_at_leaf` replaces the addressed `g[a1, ..., ak]` with
`Sequence[a1, ..., ak]`, and `flatten_sequences` (`src/eval.c`) splices that
`Sequence` into the parent on the next pass. Because each targeted slot holds a
single `Sequence` node during the walk, `arg_count` never changes mid-walk, so
several positions resolve against the *original* expression with no index-shift
bookkeeping (unlike `Insert`/`Delete`, which must sort positions descending).

**Data structures.** No bespoke structure — the walker owns the traversal and the
new tree is built by `expr_copy` of the leaf's arguments into the `Sequence`.
`flatten_at_atomic` guards the `Rational`/`Complex` nodes (stored as
`EXPR_FUNCTION` but atomic), so `FlattenAt[{1/2, x}, 1]` never manufactures
`Sequence[1, 2]`.

**Complexity / limits.** Linear in the addressed subexpressions' sizes. A visible
`NDArray` is atomic, so it is first materialised once to a nested `List`
(`ndarray_to_nested_list`) and flattened there; the result is ragged by
construction, so it is deliberately never repacked. A structural head — not on
`pack.c`'s `AWARE` list and with no `Compile[]` lowering. `Protected`.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [MapAt](../../data-structures/MapAt/), [ReplaceAt](../../assignment-and-rules/ReplaceAt/)

- Source: [`src/list/flatten_at.c`](https://github.com/stblake/mathilda/blob/main/src/list/flatten_at.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_flatten_at.c`](https://github.com/stblake/mathilda/blob/main/tests/test_flatten_at.c)

## Notes & additional examples

### Notes

`FlattenAt` removes the head of the subexpression at each position and splices
its arguments into the surrounding expression. Position resolution is the shared
`MapAt`/`ReplaceAt` walker, so `{2}` (one deep path) and `{{2}, {4}}` (two
separate paths) mean exactly what they do there, and an out-of-range position
leaves the call unevaluated. Because the head at the position is *removed*,
`FlattenAt[{1, {{2}, {3}}, 4}, 2]` gives `{1, {2}, {3}, 4}` — to flatten *within*
a part instead, use `MapAt[Flatten, ...]`. Several positions are resolved against
the original expression with no index bookkeeping, since each targeted slot
briefly holds a single `Sequence` the evaluator splices on the next pass.
