### Worked examples

```mathematica
In[1]:= Context[]
Out[1]= "Global`"

In[2]:= Context[Sin]
Out[2]= "System`"
```

### Notes

`Context[]` returns the current context (the same value as `$Context`), while `Context[sym]` returns the context in which a symbol lives — built-ins such as `Sin` reside in `"System`"`.
