### Worked examples

```mathematica
In[1]:= BooleanQ[True]
```

```mathematica
In[1]:= BooleanQ[2 > 1]  (* the inequality reduces to True first *)
```

```mathematica
In[1]:= BooleanQ[x]  (* a plain symbol is not a boolean *)
```

```mathematica
In[1]:= BooleanQ[1]  (* nor is a number *)
```

### Notes

`BooleanQ[expr]` gives `True` when `expr` is literally the symbol `True` or `False`,
and `False` for everything else. The argument is evaluated first, so `BooleanQ[2 > 1]`
sees the reduced `True`.

It tests the *symbol*, which is the key difference from `TrueQ`: `BooleanQ[False]` is
`True` (it is a boolean), whereas `TrueQ[False]` is `False` (it is not `True`).
`BooleanQ` is `Protected`.
