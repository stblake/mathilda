# Diagonal

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Diagonal[m]
    gives the list of elements on the leading diagonal of the matrix m
    (length Min[rows, cols], so it works for a non-square m).
Diagonal[m, k]
    gives the elements on the k-th diagonal of m: k > 0 above the leading
    diagonal, k < 0 below.  An out-of-range k gives {}.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Diagonal[{{a, b, c}, {d, e, f}, {g, h, i}}]
Out[1]= {a, e, i}

In[2]:= Diagonal[{{a, b, c}, {d, e, f}, {g, h, i}}, 1]
Out[2]= {b, f}

In[3]:= Diagonal[{{a, b, c}, {d, e, f}, {g, h, i}}, -1]
Out[3]= {d, h}

In[4]:= Diagonal[{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}}]
Out[4]= {1, 6, 11}

In[5]:= Diagonal[{{a, b, c, d}, {e, f, g, h}, {i, j, k, l}, {m, n, o, p}}, -2]
Out[5]= {i, n}
```

## Implementation notes

- `Protected`.
- `Diagonal[m]` gives the leading diagonal `{m[[1,1]], m[[2,2]], ...}`, of length

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
