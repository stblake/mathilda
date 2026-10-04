### Worked examples

```mathematica
In[1]:= KeyUnion[{<|a -> 1, b -> 2|>, <|b -> 3, c -> 4|>}]  (* pad both to keys a, b, c *)
```

### Notes

Every returned association carries the same keys, in first-appearance order across
the inputs, so the list is ready for tabular / row-wise processing. A key missing
from a given association is filled with `Missing["KeyAbsent", key]`; the two-argument
`KeyUnion[{…}, f]` fills it with `f[key]` instead.
