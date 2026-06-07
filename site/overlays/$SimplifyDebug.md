### Worked examples

```mathematica
In[1]:= $SimplifyDebug
Out[1]= False
```

### Notes

`$SimplifyDebug` is a global flag, default `False`. When set to `True`,
`Simplify` prints one stderr line per transform invocation
(`/Name/: <input> -> <output> [<ms> ms]`), which is useful for diagnosing slow
`Simplify` calls.
