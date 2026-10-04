### Worked examples

```mathematica
In[1]:= KeyIntersection[{<|a -> 1, b -> 2, c -> 3|>, <|b -> 20, c -> 30, d -> 40|>}]
```

### Notes

Each association is cut down to the keys shared by all of them, in the key order of
the *first*; the values come from each association unchanged, so the common keys can
differ in value across the result. `KeyComplement` is the complementary operation —
the first association's entries whose keys appear in none of the others.
