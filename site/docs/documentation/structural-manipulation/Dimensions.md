# Dimensions

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Dimensions[expr]
    gives a list of the dimensions of expr.
Dimensions[expr, n]
    gives the dimensions of expr down to at most level n.

expr is treated as a full array only at levels where every sub-piece
shares the same head and length; the walk halts at the first ragged
level. Dimensions always returns a List, including the empty List {}
for atomic expressions.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Dimensions[{{1, 2}, {3, 4}}]
Out[1]= {2, 2}

In[2]:= Dimensions[{{a, b, c}, {d, e}, {f}}]
Out[2]= {3}

In[3]:= Dimensions[{{{{a, b}}}}]
Out[3]= {1, 1, 1, 2}

In[4]:= Dimensions[{{{{a, b}}}}, 2]
Out[4]= {1, 1}

In[5]:= Dimensions[1]
Out[5]= {}
```

## Implementation notes

**Algorithm.** `builtin_dimensions` (in `src/core.c`) measures the shape of a rectangular nested structure all of whose levels share the same head (taken from the top-level expression's head). The recursive helper `get_dimensions` records each level's `arg_count`, then recurses into the first child to get the candidate sub-shape and verifies every sibling has identical depth and dimensions; as soon as the structure becomes ragged it stops and returns the dimensions found so far. An optional second argument caps the depth (`Infinity` maps to the internal cap).

**Data structures.** Fixed-size `int64_t dims[DIMENSIONS_MAX_DEPTH]` stack buffers (`DIMENSIONS_MAX_DEPTH = 64`) hold per-level extents; the result is a `List` of integers.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
