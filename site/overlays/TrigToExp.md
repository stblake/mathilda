### Worked examples

```mathematica
In[1]:= TrigToExp[Sin[x]]
Out[1]= (-1/2*I) E^(I x) + (1/2*I) E^(-I x)
```

```mathematica
In[1]:= TrigToExp[Cos[x]]
Out[1]= 1/2 E^(-I x) + 1/2 E^(I x)
```

Hyperbolic functions become real exponentials:

```mathematica
In[1]:= TrigToExp[Cosh[x]]
Out[1]= 1/2 E^x + 1/2 E^(-x)
```

The inverse functions are rewritten via their logarithmic forms — here the
classic `arctan` identity:

```mathematica
In[1]:= TrigToExp[ArcTan[x]]
Out[1]= (-1/2*I) Log[1 + I x] + (1/2*I) Log[1 - I x]
```

```mathematica
In[1]:= TrigToExp[ArcSin[x]]
Out[1]= -I Log[I x + Sqrt[1 - x^2]]
```
