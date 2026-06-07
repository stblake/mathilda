### Worked examples

```mathematica
In[1]:= a === a
Out[1]= True

In[2]:= 1 === 1.
Out[2]= False

In[3]:= f[x] === f[x]
Out[3]= True
```

### Notes

`SameQ` (`===`) is a structural test that always returns `True` or `False`. Numerically equal but distinct heads, such as `1` (Integer) and `1.` (Real), are not the same.
