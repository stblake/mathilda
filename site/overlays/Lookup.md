### Worked examples

```mathematica
In[1]:= Lookup[<|a -> 1, b -> 2|>, a]
In[2]:= Lookup[<|a -> 1, b -> 2|>, c, 0]  (* a default for the absent key *)
In[3]:= Lookup[<|a -> 1, b -> 2, c -> 3|>, {c, a}]  (* several keys at once *)
In[4]:= Lookup[<|a -> 1|>, c]  (* no default: a Missing object *)
```

### Notes

`Lookup` is `HoldAll`, so the default is evaluated only when the key is actually
absent (and once per absent key); `Lookup[<|a -> 1|>, a, Print["x"]]` prints
nothing. A list of keys shares a single index build, so it costs O(n + m) rather than
m separate scans. `Key[k]` looks up one literal key even when `k` is a list.
