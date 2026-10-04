### Worked examples

```mathematica
In[1]:= TakeLargest[{3, 1, 4, 1, 5, 9, 2, 6}, 3]  (* the three largest, descending *)
```

```mathematica
In[1]:= TakeLargest[<|a -> 3, b -> 1, c -> 5|>, 2]  (* over an association, the entries with the largest values *)
```

### Notes

`TakeLargest[list, n]` gives the `n` largest elements of `list` in descending
order, by Mathilda's canonical order. Over an association it returns the `n`
entries with the largest values, as an association. Use `TakeLargestBy` to rank by
a key function instead of by the elements themselves. If `n` exceeds the length,
all elements are returned.
