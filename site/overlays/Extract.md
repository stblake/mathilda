### Worked examples

```mathematica
In[1]:= Extract[{a,b,c},2]
Out[1]= b

In[2]:= Extract[{{a,b},{c,d}},{2,1}]
Out[2]= c

In[3]:= Extract[{{1,2},{3,4}},{{1,1},{2,2}}]
Out[3]= {1, 4}
```

### Notes

A single integer extracts one element; a position list like `{2,1}` extracts a nested part. A list of positions returns the list of extracted values.
