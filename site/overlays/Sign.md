### Worked examples

```mathematica
In[1]:= Sign[-7]
Out[1]= -1

In[2]:= Sign[0]
Out[2]= 0

In[3]:= Sign[{-2, 0, 5}]
Out[3]= {-1, 0, 1}
```

```mathematica
In[1]:= Sign[3 + 4 I]
Out[1]= 3/5 + 4/5*I
```

```mathematica
In[1]:= Sign[(1 + I)^2]
Out[1]= I

In[2]:= Sign[2 - 2 I]
Out[2]= (1/2 - 1/2*I) Sqrt[2]
```

### Notes

For real `x`, `Sign[x]` is -1, 0, or 1. For a nonzero complex `z` it returns the unit-modulus direction `z/Abs[z]`: `(1+I)^2 = 2I` points straight up the imaginary axis (`I`), while `2 - 2I` lies on the diagonal and returns the exact unit vector `(1 - I)/Sqrt[2]`. Sign is Listable.
