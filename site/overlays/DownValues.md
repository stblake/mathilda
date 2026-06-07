### Worked examples

```mathematica
In[1]:= square[n_] := n*n
Out[1]= Null

In[2]:= DownValues[square]
Out[2]= {n_^2 -> n^2}

In[3]:= square[7]
Out[3]= 49
```

### Notes

`DownValues[s]` returns the pattern rules defined on `s` via `f[args] := ...`. Each entry is the stored `lhs -> rhs` rule the evaluator tries when `s` is called.
