### Worked examples

```mathematica
In[1]:= TakeLargestBy[{1, -5, 3, -2, 4}, Abs, 2]  (* the two with the largest absolute value *)
```

```mathematica
In[1]:= TakeLargestBy[{{1, 9}, {5, 2}, {3, 7}}, Last, 2]  (* rank the pairs by their last entry *)
```

### Notes

`TakeLargestBy[list, f, n]` gives the `n` elements of `list` for which `f` is
largest, in descending order of `f`. The ranking uses `f` of each element but the
original elements are returned — here `-5` and `4` have the two largest absolute
values. Over an association it ranks by `f` of each value. If `n` exceeds the
length, all elements are returned.
