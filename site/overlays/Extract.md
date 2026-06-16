### Worked examples

```mathematica
In[1]:= Extract[{a,b,c},2]
Out[1]= b

In[2]:= Extract[{{a,b},{c,d}},{2,1}]
Out[2]= c

In[3]:= Extract[{{1,2},{3,4}},{{1,1},{2,2}}]
Out[3]= {1, 4}
```

```mathematica
In[1]:= Extract[a+b+c, 0]
Out[1]= Plus
```

```mathematica
In[1]:= Extract[x^4 + 2 x^2 + 1, Position[x^4 + 2 x^2 + 1, x]]
Out[1]= {x, x}
```

### Notes

A single integer extracts one element; a position list like `{2,1}` extracts a nested part. A list of positions returns the list of extracted values. Index `0` selects the head, so `Extract[expr, 0]` returns `Head[expr]`. Because the position syntax matches the output of `Position`, the two compose directly: feeding `Position[expr, x]` back into `Extract` pulls out every occurrence of `x`.
