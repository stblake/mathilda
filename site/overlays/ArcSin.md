### Worked examples

```mathematica
In[1]:= ArcSin[1/2]
Out[1]= 1/6 Pi

In[2]:= ArcSin[1]
Out[2]= 1/2 Pi

In[3]:= N[ArcSin[0.5]]
Out[3]= 0.523599

In[4]:= ArcSin[Sin[x]]
Out[4]= ArcSin[Sin[x]]
```

```mathematica
In[1]:= ArcSin[Sqrt[3]/2]
Out[1]= 1/3 Pi

In[2]:= ArcSin[I]
Out[2]= I ArcSinh[1]
```

```mathematica
In[1]:= D[ArcSin[x], x]
Out[1]= 1/Sqrt[1 - x^2]

In[2]:= Series[ArcSin[x], {x, 0, 7}]
Out[2]= x + 1/6 x^3 + 3/40 x^5 + 5/112 x^7 + O[x]^8
```

### Notes

`ArcSin[z]` gives the principal inverse sine, in `[-Pi/2, Pi/2]` for real `z` in `[-1, 1]`. The inverse-of-forward composition `ArcSin[Sin[x]]` is deliberately *not* folded to `x`, since that holds only on the principal branch. Exact angles such as
`ArcSin[Sqrt[3]/2] == 1/3 Pi` are recognised, and imaginary arguments map to the
hyperbolic inverse, `ArcSin[I] == I ArcSinh[1]`. Differentiation gives the
familiar `1/Sqrt[1 - x^2]`, whose Maclaurin expansion reproduces the classical
`ArcSin` series `x + x^3/6 + 3 x^5/40 + 5 x^7/112 + ...`. `ArcSin` is Listable.
