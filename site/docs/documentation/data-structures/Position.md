# Position

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Position[expr, pattern] gives a list of the positions at which objects matching pattern appear in expr.`**

**`Position[expr, pattern, levelspec] finds only objects that appear on levels specified by levelspec.`**

**`Position[expr, pattern, levelspec, n] gives the positions of the first n objects found.`**

**`Position[pattern] represents an operator form of Position that can be applied to an expression.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Position[<|"a" -> 1, "b" -> 2, "c" -> 1|>, 1]
Out[1]= {{Key["a"]}, {Key["c"]}}

In[2]:= Position[<|"a" -> {1, 2}, "b" -> 3, "c" -> 1|>, 1]
Out[2]= {{Key["a"], 1}, {Key["c"]}}

In[3]:= Extract[<|"a" -> {10, 20}|>, {Key["a"], 2}]
Out[3]= 20
```

### Applications (5)

```mathematica
In[4]:= Position[{a,b,a},a]
Out[4]= {{1}, {3}}

In[5]:= Position[{1,2,3,4},_?EvenQ]
Out[5]= {{2}, {4}}

In[6]:= Position[{1, 2, 3, 4, 5, 6}, _?PrimeQ]
Out[6]= {{2}, {3}, {5}}

In[7]:= Position[x^2 + y^2 + z^2, _Symbol]
Out[7]= {{0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {3, 0}, {3, 1}}

In[8]:= Position[Sin[Cos[x] + Tan[x]], x, Infinity]
Out[8]= {{1, 1, 1}, {1, 2, 1}}
```

## Implementation notes

**Algorithm.** `builtin_position` (`src/patterns.c`) returns the index paths of every subexpression matching the pattern. It defaults to levels `0..Infinity` and `Heads -> True` (so heads are searched and positions can contain `0`); the level-spec, `Heads` option, and optional result limit `n` are parsed as in `Cases`. The worker `do_position_at_level` carries a running `int64_t* current_path` (the accumulated index trail) and recurses depth-first pre-order: descending into the head appends index `0`, descending into argument `i` appends `i+1`. At each in-range node it calls `match(e, pattern, env)`; on success it materialises the current path as a `List` of integers and appends it to the results. Negative level-specs are resolved against `get_expr_depth_patterns`; collection stops at `max_results`. `Position[pat]` with one argument returns the operator form `Function[Position[#1, pat]]`.

**Data structures.** A reused `int64_t` path array (reallocated one deeper per recursion level), each match snapshotted into an `Expr` `List` of `EXPR_INTEGER`s; a growable `Expr**` results buffer wrapped into the final `List`; one `MatchEnv` per node.

**Attributes:** `Protected`.

## References

**See also:** [Part](../../structural-manipulation/Part/), [Extract](../../structural-manipulation/Extract/)

- Source: [`src/patterns.c`](https://github.com/stblake/mathilda/blob/main/src/patterns.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_integrate_risch_transcendental.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_risch_transcendental.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

Returns a list of position specifications (each itself a list), suitable for use with `Extract`.
