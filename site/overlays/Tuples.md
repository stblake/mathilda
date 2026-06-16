### Worked examples

```mathematica
In[1]:= Tuples[{a, b}, 2]
Out[1]= {{a, a}, {a, b}, {b, a}, {b, b}}
```

The flat-product form `Tuples[{list1, list2, ...}]` takes one element from each list — here every pairing of a coin with a die face:

```mathematica
In[1]:= Tuples[{{1, 2}, {a, b, c}}]
Out[1]= {{1, a}, {1, b}, {1, c}, {2, a}, {2, b}, {2, c}}
```

Tuples enumerates the full sample space, so it composes with `Select` for combinatorial searches. The ways two dice sum to 7:

```mathematica
In[1]:= Select[Tuples[Range[6], 2], #[[1]] + #[[2]] == 7 &]
Out[1]= {{1, 6}, {2, 5}, {3, 4}, {4, 3}, {5, 2}, {6, 1}}
```

The count of `n`-tuples from a `k`-element set is exactly `k^n`:

```mathematica
In[1]:= Length[Tuples[{a, b, c}, 4]]
Out[1]= 81
```

### Notes

`Tuples[list, n]` builds the `n`-fold Cartesian power of `list`; `Tuples[{l1, l2, ...}]` builds the heterogeneous Cartesian product, one factor per list. Results are generated in odometer (lexicographic) order, with the last index varying fastest.
