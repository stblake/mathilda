### Worked examples

```mathematica
In[1]:= TakeWhile[{2, 4, 6, 1, 8}, EvenQ]  (* stops at 1; the trailing 8 is not taken *)
In[2]:= TakeWhile[{1, 2, 3, 4}, # < 3 &]
```

### Notes

`TakeWhile` returns the leading run of passing elements and stops at the first
failure — later passing elements are not included. Over an association it tests the
values and returns an association. `LengthWhile` gives the length of the same run.
