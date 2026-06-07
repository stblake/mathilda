### Worked examples

```mathematica
In[1]:= VerifySolutions
Out[1]= VerifySolutions

In[2]:= Solve[x^2 - 3 x + 2 == 0, x, VerifySolutions -> True]
Out[2]= {{x -> 1}, {x -> 2}}
```

### Notes

`VerifySolutions` is a `Solve` option (default `Automatic`), not a function;
the bare symbol evaluates to itself. It controls whether each returned solution
is checked by back-substitution before being reported, discarding spurious
roots introduced during the solving process.
