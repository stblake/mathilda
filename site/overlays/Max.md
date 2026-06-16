### Worked examples

```mathematica
In[1]:= Max[3, 7, 2]
Out[1]= 7
```

```mathematica
In[1]:= Max[{1, 5}, {9, 2}]
Out[1]= 9
```

```mathematica
In[1]:= Max[2^100, 3^60, 5^40]
Out[1]= 1267650600228229401496703205376
```

```mathematica
In[1]:= Max[Abs[Eigenvalues[{{2, 1}, {1, 2}}]]]
Out[1]= 3
```

```mathematica
In[1]:= Max[x, 3, x]
Out[1]= Max[3, x]
```

### Notes

`Max[x1, x2, ...]` returns the numerically largest argument, and `Max` of a
collection of lists returns the largest element across all of them. Comparisons
are exact: `Max[2^100, 3^60, 5^40]` resolves a contest between three large
bignums (with automatic GMP promotion) and returns `2^100`, the actual winner.
Because `Max` is variadic and flattens lists, it composes naturally with other
operations — `Max[Abs[Eigenvalues[...]]]` computes the spectral radius of a
matrix in one line. When some arguments are non-numeric symbols, `Max` keeps the
call symbolic but still simplifies what it can: it discards duplicate operands
and drops any argument provably smaller than another, so `Max[x, 3, x]` collapses
to `Max[3, x]`.
