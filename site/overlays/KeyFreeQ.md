### Worked examples

```mathematica
In[1]:= KeyFreeQ[<|a -> 1, b -> 2|>, c]  (* true: the key c is absent *)
In[2]:= KeyFreeQ[<|a -> 1, b -> 2|>, a]  (* false: a is present *)
In[3]:= KeyFreeQ[<|1 -> x, 2 -> y|>, _Integer]  (* a pattern: some key matches, so not free *)
```

### Notes

`KeyFreeQ` is the exact complement of `KeyMemberQ`. Its second argument is a
**pattern**, so `KeyFreeQ[a, _]` is `False` for any non-empty association. This is
the one place it differs from `KeyExistsQ`, which looks for the literal key `_`.
Both an association and a bare list of rules are accepted.
