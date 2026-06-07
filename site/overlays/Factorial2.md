### Worked examples

```mathematica
In[1]:= Factorial2[7]
Out[1]= 105

In[2]:= Factorial2[8]
Out[2]= 384

In[3]:= 0!!
Out[3]= 1
```

### Notes

`Factorial2[n]` (also typeset `n!!`) is the double factorial: it multiplies down in steps of 2, so `7!! = 7*5*3*1 = 105` and `8!! = 8*6*4*2 = 384`. By convention `0!! = (-1)!! = 1`.
