### Worked examples

```mathematica
In[1]:= TakeSmallest[{3, 1, 4, 1, 5, 9, 2, 6}, 3]  (* the three smallest, ascending *)
```

```mathematica
In[1]:= TakeSmallest[<|a -> 3, b -> 1, c -> 5|>, 2]  (* over an association, the entries with the smallest values *)
```

### Notes

`TakeSmallest[list, n]` gives the `n` smallest elements of `list` in ascending
order, the mirror of `TakeLargest`. Over an association it returns the `n` entries
with the smallest values. Use `TakeSmallestBy` to rank by a key function. If `n`
exceeds the length, all elements are returned.
