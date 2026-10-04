### Worked examples

```mathematica
In[1]:= DeleteDuplicatesBy[{1, 2, 3, 4, 5, 6}, EvenQ]  (* first odd and first even survive *)
```

```mathematica
In[1]:= DeleteDuplicatesBy[{-1, 1, 2, -2, 3}, Abs]  (* keyed by |x| *)
```

### Notes

`DeleteDuplicatesBy[list, f]` keeps the first element for each distinct value of
`f[element]`, preserving order — the deduplicating cousin of `GatherBy`, returning
one representative per group instead of the groups themselves. With
`f = Abs` it collapses `±x` to whichever sign appears first. Over an association
`f` is applied to each value and the surviving entries are returned as an
association with keys preserved. Use `DeleteDuplicates` for the plain
`f = Identity` case.
