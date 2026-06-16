### Worked examples

```mathematica
In[1]:= Hypergeometric2F1[1, 1, 2, z]
Out[1]= -Log[1 - z]/z
```

A non-positive integer upper parameter terminates the Gauss series to a
polynomial — here `(1 - z)^3`:

```mathematica
In[1]:= Hypergeometric2F1[-3, 1, 1, z]
Out[1]= 1 - 3 z + 3 z^2 - z^3
```

The function reproduces the elementary inverse trig functions: `z *
2F1[1/2, 1/2, 3/2, z^2] = ArcSin[z]`. Checked at `z = 1/2` to 40 digits:

```mathematica
In[1]:= N[Hypergeometric2F1[1/2, 1/2, 3/2, 1/4]/2, 40]
Out[1]= 0.52359877559829887307710723054658381403285

In[2]:= N[ArcSin[1/2], 40]
Out[2]= 0.52359877559829887307710723054658381403285
```

### Notes

`Hypergeometric2F1[a, b, c, z]` is the Gauss hypergeometric function `HypergeometricPFQ[{a, b}, {c}, z]`, convergent for `|z| < 1` (and by termination for non-positive integer `a` or `b`). Many elementary functions are special cases: `2F1[1, 1, 2, z] = -Log[1 - z]/z` and `z * 2F1[1/2, 1/2, 3/2, z^2] = ArcSin[z]`.
