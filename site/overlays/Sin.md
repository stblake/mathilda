### Worked examples

```mathematica
In[1]:= Sin[Pi/6]
Out[1]= 1/2

In[2]:= N[Sin[1]]
Out[2]= 0.841471

In[3]:= Sin[{0, Pi/6, Pi/2}]
Out[3]= {0, 1/2, 1}

In[4]:= Sin[ArcSin[x]]
Out[4]= x
```

### Notes

The argument is in radians; rational multiples of `Pi` reduce to exact values while numeric inputs route to libm (Real) or MPFR. `Sin` is Listable, so it threads over lists element-wise.
