### Worked examples

```mathematica
In[1]:= ArcCos[0]
Out[1]= 1/2 Pi

In[2]:= ArcCos[1/2]
Out[2]= 1/3 Pi

In[3]:= ArcCos[1]
Out[3]= 0

In[4]:= N[ArcCos[0.5]]
Out[4]= 1.0472
```

```mathematica
In[1]:= ArcCos[-1/2]
Out[1]= 2/3 Pi

In[2]:= ArcCos[Sqrt[2]/2]
Out[2]= 1/4 Pi
```

```mathematica
In[1]:= N[ArcCos[2], 20]
Out[1]= 0.0 + 1.31695789692481670862*I
```

### Notes

`ArcCos[z]` gives the principal inverse cosine, in `[0, Pi]` for real `z` in `[-1, 1]`. Special exact angles are recognised: `ArcCos[-1/2]` is `2/3 Pi` and `ArcCos[Sqrt[2]/2]` is `1/4 Pi`. For `|z| > 1` the result moves off the real axis onto the branch cut, so `ArcCos[2]` is purely imaginary, `1.3169... I`. `ArcCos` is Listable.
