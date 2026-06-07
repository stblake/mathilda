### Worked examples

```mathematica
In[1]:= MatchQ[f[], f[___]]
Out[1]= True

In[2]:= MatchQ[f[1, 2], f[___]]
Out[2]= True

In[3]:= FullForm[___]
Out[3]= BlankNullSequence[]
```

### Notes

`___` (`BlankNullSequence[]`) matches a sequence of *zero or more* expressions, so it succeeds even on an empty argument list (Out[1]). Compare `__`, which requires at least one.
