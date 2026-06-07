### Worked examples

```mathematica
In[1]:= 7 - 3
Out[1]= 4

In[2]:= Subtract[7, 3]
Out[2]= 4

In[3]:= a - b - c
Out[3]= a - b - c
```

### Notes

`x - y` is rewritten to `Plus[x, Times[-1, y]]`, so subtraction inherits Plus's flattening and canonical ordering rather than existing as a distinct head in the result.
