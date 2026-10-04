### Worked examples

```mathematica
In[1]:= Keys[<|a -> 1, b -> 2, c -> 3|>]
In[2]:= Keys[<|a -> 1, b -> 2|>, f]  (* wrap each key as f[k] *)
In[3]:= Keys[{a -> 1, b -> 2}]  (* a bare list of rules is accepted too *)
```

### Notes

Keys are returned in insertion order. `Keys[assoc, f]` maps `f` over each key, and
both forms thread over a list of associations or a list of rules. The dual builtin
`Values` reads the right-hand sides instead.
