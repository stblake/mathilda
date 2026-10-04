### Worked examples

```mathematica
In[1]:= SortBy[{{2, 1}, {1, 2}, {3, 0}}, First]  (* order by the first part *)
In[2]:= SortBy[{-3, 1, -2}, Abs]  (* order by magnitude *)
In[3]:= SortBy[{<|n -> 3|>, <|n -> 1|>, <|n -> 2|>}, #n &]  (* key extracted from each record *)
```

### Notes

The key `f[element]` is computed once per element and cached, so an expensive `f`
runs `n` times, not `n log n`. Ties fall back to the canonical order of the elements
and then to original position (a stable order). A list of functions `{f1, f2, …}`
sorts by `f1`, breaking ties by `f2`, and so on; the three-argument form takes an
explicit ordering function.
