### Worked examples

```mathematica
In[1]:= Hold[1+1]
Out[1]= Hold[1 + 1]

In[2]:= ReleaseHold[Hold[1+1]]
Out[2]= 2

In[3]:= Hold[1+1, 2+2]
Out[3]= Hold[1 + 1, 2 + 2]
```

### Notes

`Hold` has attribute `HoldAll`, so it keeps every argument unevaluated; wrap an argument in `Evaluate[...]` to force a single evaluation, and use `ReleaseHold` to strip the `Hold` and evaluate the contents.
