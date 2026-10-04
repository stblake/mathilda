### Worked examples

```mathematica
In[1]:= TakeSmallestBy[{1, -5, 3, -2, 4}, Abs, 2]  (* the two with the smallest absolute value *)
```

```mathematica
In[1]:= TakeSmallestBy[{{1, 9}, {5, 2}, {3, 7}}, Last, 2]  (* rank the pairs by their last entry *)
```

### Notes

`TakeSmallestBy[list, f, n]` gives the `n` elements of `list` for which `f` is
smallest, in ascending order of `f` — the mirror of `TakeLargestBy`. The ranking
uses `f` of each element but returns the original elements. Over an association it
ranks by `f` of each value. If `n` exceeds the length, all elements are returned.
