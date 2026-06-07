### Worked examples

```mathematica
In[1]:= ByteCount[5]
Out[1]= 48

In[2]:= ByteCount[{1, 2, 3}]
Out[2]= 269

In[3]:= ByteCount[x^2 + 1]
Out[3]= 381
```

### Notes

`ByteCount` reports the number of bytes Mathilda uses internally to store an expression, including every subexpression node. The exact figures are an implementation detail of the build and are useful mainly for comparing the relative size of expressions.
