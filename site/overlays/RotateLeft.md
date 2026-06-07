### Worked examples

```mathematica
In[1]:= RotateLeft[{1,2,3,4},1]
Out[1]= {2, 3, 4, 1}

In[2]:= RotateLeft[{a,b,c,d,e},2]
Out[2]= {c, d, e, a, b}
```

### Notes

Cyclically shifts elements `n` positions to the left, wrapping the front elements to the back.
