### Worked examples

```mathematica
In[1]:= ArcTan[1]
Out[1]= 1/4 Pi

In[2]:= ArcTan[Sqrt[3]]
Out[2]= 1/3 Pi
```

```mathematica
In[1]:= ArcTan[-1, 1]
Out[1]= 3/4 Pi

In[2]:= ArcTan[{0, 1}]
Out[2]= {0, 1/4 Pi}
```

A pure-imaginary argument crosses over to the inverse hyperbolic tangent:

```mathematica
In[1]:= ArcTan[I/2]
Out[1]= I ArcTanh[1/2]
```

The Taylor series recovers the Gregory series, and the antiderivative is closed-form:

```mathematica
In[1]:= Series[ArcTan[x], {x, 0, 9}]
Out[1]= x - 1/3 x^3 + 1/5 x^5 - 1/7 x^7 + 1/9 x^9 + O[x]^10

In[2]:= Integrate[ArcTan[x], x]
Out[2]= 1/2 (2 x ArcTan[x] - Log[1 + x^2])
```

Machin's 1706 formula then evaluates Pi to 40 digits:

```mathematica
In[1]:= N[16 ArcTan[1/5] - 4 ArcTan[1/239], 40]
Out[1]= 3.1415926535897932384626433832795028841975
```

### Notes

`ArcTan[z]` gives the principal inverse tangent, in `(-Pi/2, Pi/2)`. The two-argument form `ArcTan[x, y]` is the quadrant-aware `atan2`, giving the argument of `x + I y` in `(-Pi, Pi]`. `ArcTan` is Listable.
