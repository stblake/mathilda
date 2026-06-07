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

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
