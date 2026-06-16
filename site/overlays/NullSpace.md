### Worked examples

The classic rank-deficient magic-like matrix has a one-dimensional kernel:

```mathematica
In[1]:= NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= {{1, -2, 1}}
```

A full-column-rank matrix has only the trivial null space, returned as the empty
list:

```mathematica
In[1]:= NullSpace[{{1, 0}, {0, 1}}]
Out[1]= {}
```

A wide matrix with a two-dimensional kernel — the basis vectors are scaled to
clear denominators and ordered with the rightmost free column first:

```mathematica
In[1]:= NullSpace[{{1, 2, 3, 4}, {2, 4, 6, 8}}]
Out[1]= {{-4, 0, 0, 1}, {-3, 0, 1, 0}, {-2, 1, 0, 0}}
```

NullSpace works symbolically — here the kernel is parametrized by `a`:

```mathematica
In[1]:= NullSpace[{{1, a}, {1, a}}]
Out[1]= {{-a, 1}}
```

The defining property `m . v == 0` can be checked directly:

```mathematica
In[1]:= m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}; m . First[NullSpace[m]]
Out[1]= {0, 0, 0}
```

### Notes

`NullSpace[m]` returns a basis for the null space of `m` (the vectors `v` with
`m . v == 0`). The matrix may be square or rectangular; full column rank yields
the empty list `{}`. For exact integer or rational input each basis vector is
scaled to clear denominators, so the result stays integer-valued whenever the
input is; symbolic input keeps its natural rational form. A `Method` option
selects the elimination algorithm (default `"DivisionFreeRowReduction"`, a
Bareiss-like fraction-free Gauss-Jordan).
