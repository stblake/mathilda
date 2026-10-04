### Worked examples

```mathematica
In[1]:= KeyMemberQ[<|a -> 1, b -> 2|>, a]  (* a is a key *)
In[2]:= KeyMemberQ[<|a -> 1, b -> 2|>, c]  (* c is not *)
In[3]:= KeyMemberQ[<|1 -> x, 2 -> y|>, _Integer]  (* the second argument is a pattern *)
```

### Notes

The second argument is a **pattern**: `KeyMemberQ[a, _]` is `True` for any non-empty
association, whereas `KeyExistsQ[a, _]` searches for the literal key `_`. A
pattern-free key takes the O(1) hash-index probe; a pattern falls back to a match
over each key. `KeyFreeQ` is the complement.
