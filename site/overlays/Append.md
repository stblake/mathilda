### Worked examples

```mathematica
In[1]:= Append[{1, 2, 3}, 4]
```

```mathematica
In[1]:= Append[{}, 1]  (* appending to the empty list *)
```

```mathematica
In[1]:= Append[f[a, b], c]  (* works on any head, not just List *)
```

```mathematica
In[1]:= Append[<|a -> 1, b -> 2|>, c -> 3]  (* a rule appends an entry to an association *)
```

### Notes

`Append[expr, elem]` adds `elem` as the last argument of `expr`, keeping the
original head. It is not restricted to lists: `Append[f[a, b], c]` gives
`f[a, b, c]`, and for an association a `key -> value` rule adds (or, for an
existing key, updates) an entry. The input is left unchanged — `Append` returns a
new expression — so use `AppendTo` for an in-place update of a variable.
