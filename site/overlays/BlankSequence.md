### Worked examples

```mathematica
In[1]:= MatchQ[f[1, 2], f[__]]
Out[1]= True

In[2]:= MatchQ[f[], f[__]]
Out[2]= False

In[3]:= FullForm[__]
Out[3]= BlankSequence[]
```

### Notes

`__` (`BlankSequence[]`) matches a sequence of *one or more* expressions, so it fails on the empty argument list (Out[2]). For zero-or-more use `BlankNullSequence` (`___`).
