### Worked examples

```mathematica
In[1]:= Merge[{<|a -> 1, b -> 2|>, <|a -> 10, c -> 3|>}, Total]  (* sum colliding keys *)
In[2]:= Merge[{<|a -> 1|>, <|a -> 2|>, <|b -> 3|>}, Identity]  (* keep every value in a list *)
```

### Notes

For each key, `f` receives the `List` of all values seen under it across the input
associations, in first-seen key order. Common combiners are `Total`, `Max`, `Mean`
and `Identity` (which keeps the raw value lists). The collection runs in one O(N)
hash pass.
