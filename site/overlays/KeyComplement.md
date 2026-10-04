### Worked examples

```mathematica
In[1]:= KeyComplement[{<|a -> 1, b -> 2, c -> 3|>, <|b -> 9|>}]  (* keys of the first not in the rest *)
```

### Notes

`KeyComplement[{a1, a2, ...}]` returns the entries of the first association `a1`
whose key appears in *none* of the later associations — the set difference of key
sets, carrying `a1`'s values and preserving `a1`'s order. It is the association
counterpart of `Complement`, useful for finding the records present in a baseline
but absent from every comparison set. The argument must be a non-empty list of
associations (or rule lists); the sibling `KeyIntersection` keeps the common keys
instead.
