# Position

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Position[expr, pattern] gives a list of the positions at which objects matching pattern appear in expr.
Position[expr, pattern, levelspec] finds only objects that appear on levels specified by levelspec.
Position[expr, pattern, levelspec, n] gives the positions of the first n objects found.
Position[pattern] represents an operator form of Position that can be applied to an expression.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_position` (`src/patterns.c`) returns the index paths of every subexpression matching the pattern. It defaults to levels `0..Infinity` and `Heads -> True` (so heads are searched and positions can contain `0`); the level-spec, `Heads` option, and optional result limit `n` are parsed as in `Cases`. The worker `do_position_at_level` carries a running `int64_t* current_path` (the accumulated index trail) and recurses depth-first pre-order: descending into the head appends index `0`, descending into argument `i` appends `i+1`. At each in-range node it calls `match(e, pattern, env)`; on success it materialises the current path as a `List` of integers and appends it to the results. Negative level-specs are resolved against `get_expr_depth_patterns`; collection stops at `max_results`. `Position[pat]` with one argument returns the operator form `Function[Position[#1, pat]]`.

**Data structures.** A reused `int64_t` path array (reallocated one deeper per recursion level), each match snapshotted into an `Expr` `List` of `EXPR_INTEGER`s; a growable `Expr**` results buffer wrapped into the final `List`; one `MatchEnv` per node.

- Defaults to levels `{0, Infinity}` with `Heads -> True`.
- Yields lists of indices in lexicographic order.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/patterns.c`](https://github.com/stblake/mathilda/blob/main/src/patterns.c)
- Specification: [`docs/spec/builtins/pattern-matching.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/pattern-matching.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Position[{a,b,a},a]
Out[1]= {{1}, {3}}

In[2]:= Position[{1,2,3,4},_?EvenQ]
Out[2]= {{2}, {4}}
```

```mathematica
In[1]:= Position[{1, 2, 3, 4, 5, 6}, _?PrimeQ]
Out[1]= {{2}, {3}, {5}}
```

```mathematica
In[1]:= Position[x^2 + y^2 + z^2, _Symbol]
Out[1]= {{0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {3, 0}, {3, 1}}
```

```mathematica
In[1]:= Position[Sin[Cos[x] + Tan[x]], x, Infinity]
Out[1]= {{1, 1, 1}, {1, 2, 1}}
```

### Notes

Returns a list of position specifications (each itself a list), suitable for use with `Extract`.
