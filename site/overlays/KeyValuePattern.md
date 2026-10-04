### Worked examples

```mathematica
In[1]:= MatchQ[<|a -> 1, b -> 2|>, KeyValuePattern[{a -> _}]]  (* has key a with any value *)
In[2]:= Cases[{<|x -> 1, y -> 2|>, <|x -> 5|>}, KeyValuePattern[{y -> v_}] :> v]  (* pull the y-values *)
```

### Notes

`KeyValuePattern` matches an association (or list of rules) that *contains* the given
`key -> pattern` entries, in any order and ignoring extra keys. Value patterns may
bind, so `KeyValuePattern[{"a" -> v_}]` captures the value at `"a"`. Requirements that
share a bound variable are resolved consistently by backtracking.
