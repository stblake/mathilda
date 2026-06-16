### Worked examples

```mathematica
In[1]:= MatchQ[{1, 1, 1}, {Repeated[1]}]
Out[1]= True

In[2]:= MatchQ[{1, 1, 2}, {Repeated[1]}]
Out[2]= False

In[3]:= MatchQ[{1, 2, 3}, {_Integer ..}]
Out[3]= True
```

The bound can be a range `{min, max}` or an exact count `{n}`. Combined with
`Cases`, this filters a list of lists by length — keeping only those with two or
three matching elements:

```mathematica
In[1]:= Cases[{{1}, {1, 1}, {1, 1, 1}, {1, 1, 1, 1}}, {Repeated[1, {2, 3}]}]
Out[1]= {{1, 1}, {1, 1, 1}}

In[2]:= MatchQ[{1, 2, 3, 4}, {Repeated[_Integer, {4}]}]
Out[2]= True
```

### Notes

`Repeated[p]` (postfix `p..`) matches a sequence of *one or more* expressions each satisfying `p`. The optional second argument bounds the count, e.g. `Repeated[p, {2, 4}]` for between two and four.
