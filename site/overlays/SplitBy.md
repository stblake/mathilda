### Worked examples

```mathematica
In[1]:= SplitBy[{1, 3, 5, 2, 4, 7}, EvenQ]  (* runs of equal parity, adjacency only *)
```

```mathematica
In[1]:= SplitBy[{1, 2, 4, 3, 5}, EvenQ]  (* a new run starts each time the key changes *)
```

### Notes

`SplitBy[list, f]` splits `list` into runs of consecutive elements that share the same value
of `f[element]`. Only **adjacent** elements are grouped: in the first example the leading
odds `{1, 3, 5}` form one run, the evens `{2, 4}` the next, and the trailing odd `7` its own
run — unlike `GatherBy`, which would collect all the odds together regardless of position.

It differs from `Split` in comparing the evaluated keys `f[e]` rather than the elements
themselves. `SplitBy[list, {f1, f2, ...}]` splits by `f1`, then splits each run by `f2`, and
so on, nesting one level deeper per function. Element order is preserved within every run.
