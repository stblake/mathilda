### Worked examples

```mathematica
In[1]:= GroupBy[{1, 2, 3, 4, 5, 6}, EvenQ]
```

```mathematica
In[1]:= GroupBy[Range[6], Mod[#, 3] &]
```

```mathematica
In[1]:= GroupBy[{1, 2, 3, 4, 5, 6}, EvenQ, Total]  (* a reducer summarises each group *)
```

```mathematica
In[1]:= GroupBy[{"apple", "pear", "plum", "fig"}, StringLength]
```

### Notes

`GroupBy[list, f]` returns `<|f[x] -> {elements}, ...|>`, grouping the elements by
the value of the key function and preserving first-key order. The three-argument
form `GroupBy[list, f, red]` applies a reducer to each group — `GroupBy[data,
key, Total]` is a group-and-sum in one step. A `keyfn -> valfn` second argument
groups by one function but collects another, and a list of key functions groups
into nested associations. It is the keyed counterpart of `GatherBy`; over an
association the entries are grouped by `f[value]` with keys preserved.
