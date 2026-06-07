### Worked examples

```mathematica
In[1]:= Cos[Pi/3]
Out[1]= 1/2

In[2]:= N[Cos[2]]
Out[2]= -0.416147

In[3]:= Cos[{0, Pi/3, Pi}]
Out[3]= {1, 1/2, -1}

In[4]:= Cos[ArcCos[x]]
Out[4]= x
```

### Notes

The argument is in radians; rational multiples of `Pi` reduce to exact values while numeric inputs route to libm / MPFR. `Cos` is Listable.
