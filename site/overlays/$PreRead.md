### Worked examples

```mathematica
In[1]:= $PreRead = (StringJoin["(", #, ")^2"] &)
Out[1]= StringJoin["(", #1, ")^2"] &

In[2]:= 3 + 4
Out[2]= 49
```

### Notes

`$PreRead`, if set, is applied to the raw text of every input expression before
it is parsed; its value must be a function taking one string and returning a
string. Here it wraps each input as `(...)^2`, so the text `3 + 4` is read as
`(3 + 4)^2`. A non-string return falls back to the original input with a
`$PreRead::strret` diagnostic.
