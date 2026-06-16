### Worked examples

Split a list into non-overlapping blocks of a given length (trailing
under-filled elements are dropped):

```mathematica
In[1]:= Partition[{a, b, c, d, e, f}, 2]
Out[1]= {{a, b}, {c, d}, {e, f}}

In[2]:= Partition[{1, 2, 3, 4, 5}, 2]
Out[2]= {{1, 2}, {3, 4}}
```

The offset `d` turns `Partition` into a sliding window — `d = 1` gives every
consecutive overlapping pair:

```mathematica
In[1]:= Partition[{1, 2, 3, 4, 5}, 2, 1]
Out[1]= {{1, 2}, {2, 3}, {3, 4}, {4, 5}}
```

Combined with `Map`, non-overlapping blocks let you reduce a stream in chunks —
the block sums of `1..12` taken four at a time:

```mathematica
In[1]:= Map[Total, Partition[Range[12], 4]]
Out[1]= {10, 26, 42}
```

A length-2 sliding window computes first differences. Applied to the squares
`1, 4, 9, 16, 25` it recovers the classical identity that consecutive squares
differ by successive odd numbers:

```mathematica
In[1]:= Map[(#[[2]] - #[[1]] &), Partition[{1, 4, 9, 16, 25}, 2, 1]]
Out[1]= {3, 5, 7, 9}
```

### Notes

`Partition[list, n]` cuts `list` into consecutive non-overlapping length-`n`
sublists, discarding a trailing remainder that cannot fill a full block.
`Partition[list, n, d]` advances by offset `d` between successive sublists:
`d = n` reproduces the non-overlapping blocks, while smaller `d` produces
overlapping moving windows (`d = 1` slides one element at a time). Pairing
`Partition` with `Map`/`Total` is the idiomatic way to express block reductions
and finite-difference / sliding-window computations.
