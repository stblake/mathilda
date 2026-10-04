### Worked examples

```mathematica
In[1]:= Values[<|a -> 1, b -> 2, c -> 3|>]
In[2]:= Values[<|a -> 1, b -> 2|>, f]  (* wrap each value as f[v] *)
```

### Notes

Values are returned in the association's insertion order, so `Keys` and `Values`
line up positionally. `Values[assoc, f]` maps `f` over each value, and both forms
thread over a list of associations or rules.
