### Worked examples

```mathematica
In[1]:= KeyExistsQ[<|a -> 1, b -> 2|>, a]
```

```mathematica
In[1]:= KeyExistsQ[<|a -> 1, b -> 2|>, c]
```

### Notes

`KeyExistsQ[assoc, key]` returns `True` if the association has the given key, else
`False`. The key is taken *literally*, not as a pattern — this is what
distinguishes `KeyExistsQ` from `KeyMemberQ`/`KeyFreeQ`, which match their second
argument as a pattern (`KeyMemberQ[a, _]` is `True` for any non-empty `a`, while
`KeyExistsQ[a, _]` looks for the literal key `_`). The test is a single `O(1)`
amortised probe through the association's key index, and it also accepts a bare
list of rules.
