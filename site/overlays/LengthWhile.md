### Worked examples

```mathematica
In[1]:= LengthWhile[{1, 2, 3, -1, 5}, Positive]  (* leading run stops at -1 *)
In[2]:= LengthWhile[{2, 4, 6, 7, 8}, EvenQ]  (* stops at 7, even though 8 follows *)
```

### Notes

Only the *leading* run counts: the scan halts at the first element failing the
predicate, so later passing elements are not counted. A numeric predicate takes a
compiled buffer scan that builds nothing. `TakeWhile` returns that leading run as a
collection instead of its length.
