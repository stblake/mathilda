### Worked examples

```mathematica
In[1]:= Pick[{a, b, c, d}, {True, False, True, False}]  (* two-argument form: keep where True *)
In[2]:= Pick[{1, 2, 3, 4, 5}, {1, 0, 1, 0, 1}, 1]  (* keep where the selector matches 1 *)
In[3]:= Pick[{1, 2, 3, 4, 5, 6}, {1, 2, 3, 4, 5, 6}, _?EvenQ]  (* selector matched against a pattern *)
```

### Notes

The selector must mirror the structure of the first argument; `Pick` walks both in
step and keeps each element whose selector matches the pattern (`True` by default).
The result keeps the first argument's head. A shape mismatch leaves `Pick[…]`
unevaluated rather than erroring.
