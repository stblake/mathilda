### Worked examples

```mathematica
In[1]:= KeyMap[f, <|a -> 1, b -> 2|>]  (* transform each key, values kept *)
In[2]:= KeyMap[EvenQ, <|1 -> a, 2 -> b, 3 -> c, 4 -> d|>]  (* distinct keys may collapse; last value wins *)
```

### Notes

Only the keys are transformed; the values ride along. If `f` sends two keys to the
same result the association re-canonicalises and the later entry's value survives, so
the output can be shorter than the input (as in the `EvenQ` example, where four keys
collapse to `False` and `True`).
