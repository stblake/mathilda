### Worked examples

```mathematica
In[1]:= GatherBy[{1, 2, 3, 4, 5, 6}, EvenQ]
```

```mathematica
In[1]:= GatherBy[Range[9], Mod[#, 3] &]
```

```mathematica
In[1]:= GatherBy[{{1, a}, {2, b}, {1, c}}, First]  (* group rows by their first column *)
```

### Notes

`GatherBy[list, f]` gathers elements with equal `f[x]` into sublists, in
first-appearance order — like `GroupBy`, but returning the groups as a plain list
of lists with the group keys dropped. It is the natural tool for grouping records
by a field (`GatherBy[rows, First]`) or partitioning numbers by a residue class.
Over an association the entries are gathered by `f[value]` into sub-associations.
Use `GroupBy` when you want the keys kept, or `Counts`/`CountsBy` when you only
need the sizes.
