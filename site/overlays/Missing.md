### Worked examples

```mathematica
In[1]:= Missing["KeyAbsent", x]  (* an inert marker: stays as written *)
In[2]:= Lookup[<|a -> 1|>, b]  (* an absent key produces one *)
In[3]:= MissingQ[Missing["Unknown"]]  (* test with MissingQ *)
```

### Notes

`Missing[…]` is a placeholder for data that is not there, not something that
evaluates. Absent-key lookups yield `Missing["KeyAbsent", key]`, and `KeyUnion` /
`JoinAcross` fill gaps with `Missing[…]`. Detect it with `MissingQ`, strip it with
`DeleteMissing`, or avoid it entirely by giving `Lookup` a default.
