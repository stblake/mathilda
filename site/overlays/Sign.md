### Worked examples

```mathematica
In[1]:= Sign[-7]
Out[1]= -1

In[2]:= Sign[0]
Out[2]= 0

In[3]:= Sign[3 + 4 I]
Out[3]= 3/5 + 4/5*I

In[4]:= Sign[{-2, 0, 5}]
Out[4]= {-1, 0, 1}
```

### Notes

For real `x`, `Sign[x]` is -1, 0, or 1. For a nonzero complex `z` it returns the unit-modulus direction `z/Abs[z]`. Sign is Listable.
