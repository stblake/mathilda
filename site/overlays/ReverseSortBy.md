### Worked examples

```mathematica
In[1]:= ReverseSortBy[{1, -5, 3, -2}, Abs]  (* descending by absolute value *)
```

```mathematica
In[1]:= ReverseSortBy[{{1, 2}, {3, 1}, {2, 5}}, Last]  (* order the pairs by their last entry, descending *)
```

### Notes

`ReverseSortBy[list, f]` sorts by `f` in descending order — the `By` companion to
`ReverseSort`, and the descending counterpart of `SortBy`. The key function `f` is
applied to each element (or each value, for an association) and the elements are
ordered by the canonical order of those keys, largest first. A three-argument
`ReverseSortBy[list, f, p]` ranks the keys with the reversed ordering function `p`.
