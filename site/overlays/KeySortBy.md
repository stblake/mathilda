### Worked examples

```mathematica
In[1]:= KeySortBy[<|"bbb" -> 1, "a" -> 2, "cc" -> 3|>, StringLength]  (* order keys by length *)
```

### Notes

Entries are ordered by `f` applied to each key, with ties keeping the association's
original order (the sort is stable). `f` is evaluated once per key. Contrast
`KeySort` (sort by the keys directly) and `SortBy` (sort by a function of the
values).
