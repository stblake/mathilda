### Worked examples

```mathematica
In[1]:= SetAttributes[r, OneIdentity]
Out[1]= Null

In[2]:= Attributes[r]
Out[2]= {OneIdentity}

In[3]:= MatchQ[a, r[a]]
Out[3]= True
```

### Notes

`OneIdentity` makes `r[x]`, `r[r[x]]`, etc. equivalent to `x` for pattern matching, so a bare `a` matches the pattern `r[a]`. It does not change evaluation — only how patterns unify against the head.
