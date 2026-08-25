# Level

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Level[expr, levelspec]`**

gives a list of all subexpressions of expr on levels specified by levelspec.

**`Level[expr, levelspec, f]`**

applies f to the sequence of subexpressions.

**`Level[expr, {-1}] gives a list of all "atomic" objects in expr.`**

<details>
<summary>Notes</summary>

Level uses standard level specifications: n            levels 1 through n Infinity     levels 1 through Infinity {n}          level n only {n1, n2}     levels n1 through n2 A positive level n consists of all parts of expr specified by n indices. A negative level -n consists of all parts of expr with depth n. Level 0 corresponds to the whole expression. With the option setting Heads-\>True, Level includes heads of expressions and their parts. Level traverses expressions in depth-first order, so that the subexpressions in the final list are ordered lexicographically by their indices.

</details>

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Level[a + f[x, y^n], {-1}]
Out[1]= {a, x, y, n}

In[2]:= Level[a + f[x, y^n], 2]
Out[2]= {a, x, y^n, f[x, y^n]}
```

### Options (1)

```mathematica
In[3]:= Level[x^2 + y^3, 3, Heads -> True]
Out[3]= {Plus, Power, x, 2, x^2, Power, y, 3, y^3}
```

### Applications (6)

```mathematica
In[4]:= Level[a + b c, {1}]
Out[4]= {a, b c}

In[5]:= Level[a + b c, {2}]
Out[5]= {b, c}

In[6]:= Level[f[g[h[x]]], {-1}]
Out[6]= {x}

In[7]:= Level[(1 + x)^2 + y, {-2}]
Out[7]= {1 + x}

In[8]:= Level[{{a, b}, {c, {d, e}}}, Infinity]
Out[8]= {a, b, {a, b}, c, d, e, {d, e}, {c, {d, e}}}

In[9]:= Level[f[g[x], y], 2, Heads -> True]
Out[9]= {f, g, x, g[x], y}
```

## Implementation notes

**Algorithm.** `builtin_level` collects the subexpressions of `e` that lie at the depths
selected by the level spec (an integer `n` = levels `1..n`, `Infinity`, `{n}` = exactly level
`n`, or `{min, max}`). It descends with the recursive `level_rec`, tracking the current depth
and appending matching nodes into a growable `Expr**` buffer, then wraps them in a `List`. The
option `Heads -> True` additionally includes function heads as level elements.

- `Protected`.
- Default option: `Heads -> False`.
- Standard level specifications:
  - `n`: levels 1 through `n`.
  - `Infinity`: levels 1 through `Infinity`.
  - `{n}`: level `n` only.
  - `{n1, n2}`: levels `n1` through `n2`.
- Positive level `n` refers to distance from the top (level 0 is the whole expression).
- Negative level `-n` refers to distance from the bottom (depth `n`).
- Level `-1` corresponds to atomic objects.
- Lists subexpressions in post-order (depth-first), resulting in lexicographic ordering of indices.

**Attributes:** `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_poly.c`](https://github.com/stblake/mathilda/blob/main/tests/test_poly.c)
- Tests: [`tests/test_solve_integers.c`](https://github.com/stblake/mathilda/blob/main/tests/test_solve_integers.c)

## Notes & additional examples

### Notes

`Level[expr, levelspec]` returns the subexpressions of `expr` on the specified
levels. A positive level `n` selects parts reachable by `n` indices; a negative
level `-n` selects parts of depth `n`; `{-1}` gives all atoms and level `0` is
the whole expression. `Level[expr, levelspec, f]` applies `f` to the sequence of
subexpressions, and `Heads -> True` includes expression heads in the traversal.
