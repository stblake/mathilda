### Worked examples

```mathematica
In[1]:= ArcTan[1]
Out[1]= 1/4 Pi

In[2]:= N[ArcTan[1]]
Out[2]= 0.785398

In[3]:= ArcTan[-1, 1]
Out[3]= 3/4 Pi

In[4]:= ArcTan[{0, 1}]
Out[4]= {0, 1/4 Pi}
```

### Notes

`ArcTan[z]` gives the principal inverse tangent, in `(-Pi/2, Pi/2)`. The two-argument form `ArcTan[x, y]` is the quadrant-aware `atan2`, giving the argument of `x + I y` in `(-Pi, Pi]`. `ArcTan` is Listable.
