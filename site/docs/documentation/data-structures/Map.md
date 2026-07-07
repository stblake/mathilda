# Map

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
f /@ expr or Map[f, expr]
    applies f to each element at level 1 of expr, preserving expr's head.
Map[f, expr, levelspec]
    applies f at the parts of expr selected by levelspec (e.g. {2} for
    level 2 only, Infinity for every level).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Map[#^2 &, <|"x" -> 3, "y" -> 4|>]
Out[1]= <|"x" -> 9, "y" -> 16|>

In[2]:= Select[<|"a" -> 1, "b" -> 2, "c" -> 3|>, # > 1 &]
Out[2]= <|"b" -> 2, "c" -> 3|>
```

## Implementation notes

**Algorithm.** `builtin_map` applies `f` to subexpressions of `expr` at the
levels selected by an optional level-spec (default `{1,1}`, the immediate
arguments). The recursion `map_at_level` works **bottom-up**: for an
`EXPR_FUNCTION` it first rebuilds the node by mapping into each argument (and the
head too when `Heads -> True`), then — if the node's current level is within
`[spec.min, spec.max]` (negative levels measured against `get_depth`) — wraps it
in `f[...]` and calls `evaluate()`. Atoms are copied, and tested for membership
of the level range only by their depth.

**Level / option parsing.** `parse_level_spec` reads an integer `n`, `{n}`,
`{m,n}`, or `Infinity`; `parse_options` reads a trailing `Heads -> True`. A
`Rule`-headed third argument is treated as an option rather than a level-spec.

**Data structures.** Pure `Expr`-tree traversal; new nodes built with
`expr_new_function`. Map, MapAll, and MapAt all share this module and the
`LevelSpec { min, max, heads }` struct.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Harold Abelson and Gerald Jay Sussman, *Structure and Interpretation of Computer Programs*, 2nd ed., §2.2.1 (sequence mapping).
- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Map[f, {a, b, c}]
Out[1]= {f[a], f[b], f[c]}
```

```mathematica
In[1]:= #^2 & /@ {1, 2, 3, 4}
Out[1]= {1, 4, 9, 16}
```

```mathematica
In[1]:= Map[Reverse, {{1, 2}, {3, 4}}]
Out[1]= {{2, 1}, {4, 3}}
```

```mathematica
In[1]:= Map[f, {{a}, {b}}, {2}]
Out[1]= {{f[a]}, {f[b]}}
```

```mathematica
In[1]:= Map[Total, {{1, 2, 3}, {4, 5, 6}}]
Out[1]= {6, 15}
```

```mathematica
In[1]:= Map[#^2 &, x + y + z]
Out[1]= x^2 + y^2 + z^2
```

### Notes

`f /@ expr` is the operator shorthand for `Map[f, expr]`. By default the function
is applied at level 1, i.e. to the immediate elements; a level specification such
as `{2}` reaches deeper into nested lists. `Map` works on any expression, not
only `List` — the head is preserved while each argument is wrapped by `f`. Pure
functions (`#^2 &`) are the idiomatic first argument.
