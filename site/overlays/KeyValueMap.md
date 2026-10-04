### Worked examples

```mathematica
In[1]:= KeyValueMap[f, <|a -> 1, b -> 2|>]  (* f sees each key and value *)
In[2]:= KeyValueMap[#1 -> #2^2 &, <|a -> 2, b -> 3|>]  (* rebuild rules from key and value *)
```

### Notes

`f` is applied to each key and value together, and the results form a plain `List`
(not an association). This is the two-argument counterpart of mapping over an
association, which sees only the values.
