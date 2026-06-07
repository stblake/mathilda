### Worked examples

```mathematica
In[1]:= MatchQ[{1, 1, 1}, {Repeated[1]}]
Out[1]= True

In[2]:= MatchQ[{1, 1, 2}, {Repeated[1]}]
Out[2]= False

In[3]:= MatchQ[{1, 2, 3}, {_Integer ..}]
Out[3]= True
```

### Notes

`Repeated[p]` (postfix `p..`) matches a sequence of *one or more* expressions each satisfying `p`. The optional second argument bounds the count, e.g. `Repeated[p, {2, 4}]` for between two and four.
