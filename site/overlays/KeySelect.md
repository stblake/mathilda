### Worked examples

```mathematica
In[1]:= KeySelect[<|1 -> a, 2 -> b, 3 -> c, 4 -> d|>, EvenQ]  (* keep even-numbered keys *)
In[2]:= KeySelect[<|1 -> a, 2 -> b, 3 -> c, 4 -> d|>, OddQ]
```

### Notes

The predicate is applied to each key; an entry is kept only when the result is
exactly `True`. Order is preserved. To filter by the values rather than the keys, use
`Select` over the association.
