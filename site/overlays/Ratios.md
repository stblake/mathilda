### Worked examples

```mathematica
In[1]:= Ratios[{1, 2, 4, 8, 16}]
Out[1]= {2, 2, 2, 2}
```

A constant ratio list is the fingerprint of a geometric sequence. Applied to the
Fibonacci numbers, `Ratios` produces the classic convergents of the golden
ratio, which numerically close in on `φ`:

```mathematica
In[1]:= Ratios[{1, 1, 2, 3, 5, 8, 13, 21}]
Out[1]= {1, 2, 3/2, 5/3, 8/5, 13/8, 21/13}

In[2]:= N[Ratios[{1, 1, 2, 3, 5, 8, 13, 21, 34, 55}]]
Out[2]= {1.0, 2.0, 1.5, 1.66667, 1.6, 1.625, 1.61538, 1.61905, 1.61765}
```

The successive ratios of factorials collapse to the integers, and `FoldList`
with `Times` reconstructs the original list from a seed and its ratios:

```mathematica
In[1]:= Ratios[{1, 1, 2, 6, 24, 120}]
Out[1]= {1, 2, 3, 4, 5}

In[2]:= FoldList[Times, 3, Ratios[{3, 6, 18, 36}]]
Out[2]= {3, 6, 18, 36}
```

### Notes

`Ratios[list]` gives the successive ratios `list[[k+1]]/list[[k]]` (length
`l - 1`). `Ratios[list, n]` iterates `n` times; `Ratios[list, {n1, n2, ...}]`
takes the `n_k`-th ratios at level `k` of a nested list. `FoldList[Times, x,
Ratios[list]]` inverts `Ratios`.
