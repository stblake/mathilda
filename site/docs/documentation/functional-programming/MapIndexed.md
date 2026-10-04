# MapIndexed

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MapIndexed[f, expr]`**

Applies f to the elements of expr, giving the part specification of each element as a second argument to f: {f\[e1, {1}\], f\[e2, {2}\], ...}.

**`MapIndexed[f, expr, levelspec]`**

Applies f to all parts of expr on the levels specified by levelspec: n            levels 1 through n Infinity     levels 1 through Infinity {n}          level n only {n1, n2}     levels n1 through n2 The default is {1}. A positive level n consists of all parts specified by n indices; a negative level -n consists of all parts with depth n, so level -1 is the atoms. Level 0 is the whole expression, whose position is {}. The position handed to f is the one Part and Extract take, so Extract\[expr, #2\] is #1. Over an association the position of a value is {Key\[k\]}, and keys are preserved. With Heads -\> True the function is applied to heads as well, a head having index 0 in its position. MapIndexed always effectively constructs a complete new expression and then evaluates it.

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= MapIndexed[f, {10, 20, 30}]
Out[1]= {f[10, {1}], f[20, {2}], f[30, {3}]}

In[2]:= MapIndexed[First[#2] + f[#1] &, {a, b, c, d}]
Out[2]= {1 + f[a], 2 + f[b], 3 + f[c], 4 + f[d]}

In[3]:= MapIndexed[f, {{a, b}, {c, d, e}}, {2}]
Out[3]= {{f[a, {1, 1}], f[b, {1, 2}]}, {f[c, {2, 1}], f[d, {2, 2}], f[e, {2, 3}]}}

In[4]:= MapIndexed[f, {{a, b}, {c, d, {e}}}, {-1}]
Out[4]= {{f[a, {1, 1}], f[b, {1, 2}]}, {f[c, {2, 1}], f[d, {2, 2}], {f[e, {2, 3, 1}]}}}

In[5]:= MapIndexed[f, h0[h1[h2[h3[h4[a]]]]], {2, -3}]
Out[5]= h0[h1[f[h2[f[h3[h4[a]], {1, 1, 1}]], {1, 1}]]]

In[6]:= MapIndexed[f, {a, b}, {0, 1}]
Out[6]= f[{f[a, {1}], f[b, {2}]}, {}]

In[7]:= MapIndexed[f, <|"a" -> 10, "b" -> 20|>]
Out[7]= <|"a" -> f[10, {Key["a"]}], "b" -> f[20, {Key["b"]}]|>

In[8]:= MapIndexed[h, <|a -> <|b -> c, p -> <|q -> r|>|>, d -> {e}|>, {2}]
Out[8]= <|a -> <|b -> h[c, {Key[a], Key[b]}], p -> h[<|q -> r|>, {Key[a], Key[p]}]|>, d -> {h[e, {Key[d], 1}]}|>
```

### Options (1)

```mathematica
In[9]:= MapIndexed[f, p[x][a, b, c], Infinity, Heads -> True]
Out[9]= f[f[p, {0, 0}][f[x, {0, 1}]], {0}][f[a, {1}], f[b, {2}], f[c, {3}]]
```

### Applications (4)

The index arrives as a LIST, {1}, {2}, ...

```mathematica
In[10]:= MapIndexed[f, {a, b, c}]
Out[10]= {f[a, {1}], f[b, {2}], f[c, {3}]}
```

Keep only the position of each element

```mathematica
In[11]:= MapIndexed[#2 &, {a, b, c}]
Out[11]= {{1}, {2}, {3}}
```

Pair each element with its integer index

```mathematica
In[12]:= MapIndexed[{#1, First[#2]} &, {x, y, z}]
Out[12]= {{x, 1}, {y, 2}, {z, 3}}
```

At level 2 the position is a two-index list

```mathematica
In[13]:= MapIndexed[f, {{a, b}, {c, d}}, {2}]
Out[13]= {{f[a, {1, 1}], f[b, {1, 2}]}, {f[c, {2, 1}], f[d, {2, 2}]}}
```

## Implementation notes

**Algorithm.** `builtin_mapindexed` wraps each selected part of `expr` in
`f[part, {position}]`, where `position` is the list of indices locating the part.
It parses an optional level spec (default `{1}`) with `parse_level_spec_strict`
and a trailing `Heads -> True` option, then rebuilds the tree bottom-up in
`mi_at_level`: each child is visited at `level + 1` carrying an `MIPath` frame
that records its index component, and a node at a level inside the spec is wrapped
with its accumulated position (`mi_position` reverses the linked frames into a
`List`). An association's parts are its *values*, positioned by `Key[k]`, and its
head is never mapped. `Rational` and `Complex` are atomic here, as they are for
`Depth` and `Level`.

A visible `NDArray` is atomic to the generic traversal, so the default level `{1}`
iterates its leading axis directly (`mi_ndarray_axis`); any other spec materialises
the array to a nested list first and repacks nothing. An empty level range (e.g.
`{3, 1}`) selects nothing and returns the expression itself, which keeps a packed
array packed.

**Data structures.** Positions are built from a stack-allocated `MIPath` linked
list (`{path, component, length}` frames) threaded through the recursion, so no
index vector is copied per node. The original subtree depth is accumulated into
`*out_depth` as the recursion unwinds — needed only by a negative level bound —
keeping the whole traversal `O(n)` rather than recomputing a depth at each node.

**Complexity / limits.** `O(n)` over the parts at a non-negative spec, with early
pass-through once past the maximum level so the common `MapIndexed[f, list]` never
walks each element's whole subtree. Negative-bound specs keep the full depth-aware
descent.

**Attributes:** `Protected`.

## References

**See also:** [Part](../../data-structures/Part/), [Extract](../../structural-manipulation/Extract/), [Rule](../../assignment-and-rules/Rule/), [RuleDelayed](../../assignment-and-rules/RuleDelayed/), [Association](../../data-structures/Association/), [List](../../other-advanced/List/)

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_assoc_atomicity.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_atomicity.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_core_algebra.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core_algebra.c)
- Tests: [`tests/test_map_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_map_ndarray.c)

## Notes & additional examples

### Notes

`MapIndexed[f, expr]` is `Map` that also hands `f` the position of each element as
a second argument: `{f[e1, {1}], f[e2, {2}], ...}`. The position is always a list
of indices, so at deeper levels it has one entry per level (`{2}` makes the
positions `{i, j}`). Inside a pure function `#1` is the element and `#2` its
position — `First[#2]` recovers the plain integer index.

A level spec selects which parts are wrapped (default `{1}`), and `Heads -> True`
also indexes heads. Over an association the values are mapped and positioned by
`Key[k]`.
