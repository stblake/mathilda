### Worked examples

```mathematica
In[1]:= StringTake["abcdef", 3]  (* the first three characters *)
```

```mathematica
In[1]:= StringTake["abcdef", -2]  (* the last two *)
```

```mathematica
In[1]:= StringTake["abcdef", {2, 4}]  (* characters 2 through 4 *)
```

### Notes

`StringTake` slices a string by byte: a positive `n` takes the first `n`, a
negative `-n` the last `n`, `{m, n}` the inclusive range, and `{m, n, s}` a
stepped range. `UpTo[n]` clamps to the available length.

Negative endpoints count from the end; indexing is 1-based. A first argument that
is a list of strings is handled per element. An out-of-range or non-integer spec
leaves the call unevaluated.
