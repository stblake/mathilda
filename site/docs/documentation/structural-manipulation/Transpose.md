# Transpose

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Transpose[list]
    Transposes the first two levels of list (swaps rows and columns of a matrix).
Transpose[list, {n1, n2, ...}]
    Gives the transpose of list so that level k in list is level nk in the result.
    The spec must be a permutation of {1, ..., r} where r is the depth of list.
    A repeated index (e.g. {1, 1}) selects the corresponding diagonal.
    list must be a rectangular array.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Transpose[{{a, b}, {c, d}}]
Out[1]= {{a, c}, {b, d}}

In[2]:= Transpose[{{a, b}, {c, d}}, {1, 1}]
Out[2]= {a, d}
```

## Implementation notes

- `Protected`.
- Works only on rectangular arrays.
- `Transpose[m, {1, 1}]` extracts the diagonal of a square matrix.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed., Cambridge University Press, 2013 — the matrix transpose and index permutations of tensors.
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Transpose[{{1, 2, 3}, {4, 5, 6}}]
Out[1]= {{1, 4}, {2, 5}, {3, 6}}
```

```mathematica
In[1]:= Transpose[{{1, 2}, {3, 4}}]
Out[1]= {{1, 3}, {2, 4}}
```

```mathematica
In[1]:= Transpose[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {1, 1}]
Out[1]= {1, 5, 9}
```

### Notes

With one argument `Transpose` swaps the first two levels of a list, turning a `2x3` matrix into a `3x2` one. The optional permutation spec generalises this to arbitrary index reorderings of a rectangular array. A repeated index in the spec — `{1, 1}` in the third example — extracts the corresponding diagonal, here the main diagonal `{1, 5, 9}` of the `3x3` matrix. The spec must be a permutation of `{1, ..., r}` where `r` is the depth of the list.
