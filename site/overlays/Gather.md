### Worked examples

```mathematica
In[1]:= Gather[{1, 7, 3, 7, 2, 3, 9}]
```

```mathematica
In[1]:= Gather[{a, b, a}]
```

```mathematica
In[1]:= Gather[{3, 1, 3, 2, 1}]
```

### Notes

`Gather[list]` partitions the elements into sublists of identical elements: two
elements share a sublist when they are structurally equal, sublists appear in
order of first occurrence, and input order is kept within each sublist. Unlike
`Split`, the grouping is not limited to adjacent runs — equal elements anywhere in
the list land together. It is the identity case of `GatherBy` (`Gather[l]` is
`GatherBy[l, Identity]`) and is itself the key-less cousin of `Tally`, which
reports counts rather than the grouped elements.
