### Worked examples

```mathematica
In[1]:= KeyTake[<|a -> 1, b -> 2, c -> 3|>, {a, c}]  (* keep only a and c *)
```

### Notes

The kept entries stay in the association's own order, not the order of the key list.
`KeyDrop` is the complement (remove the listed keys). Both thread over a list of
associations, which is the column-of-records form.
