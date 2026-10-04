### Worked examples

```mathematica
In[1]:= KeyDrop[<|a -> 1, b -> 2, c -> 3|>, b]
```

```mathematica
In[1]:= KeyDrop[<|a -> 1, b -> 2, c -> 3|>, {a, c}]  (* drop several keys at once *)
```

```mathematica
In[1]:= KeyDrop[{<|a -> 1, b -> 2|>, <|a -> 3, b -> 4|>}, a]  (* threads over a list of records *)
```

### Notes

`KeyDrop[assoc, key]` returns the association with the given key removed, keeping
the order of the remaining entries; the second argument may be a single key or a
list of keys. Given a list of associations it threads — dropping the key(s) from
each record, which is the column-of-records form. It is the complement of
`KeyTake`, and both build a hash index of the key set once so the work is linear
rather than quadratic in the key count.
