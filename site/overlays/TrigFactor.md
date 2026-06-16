### Worked examples

`TrigFactor` collapses an expanded angle-addition expression back to a single
trigonometric function:

```mathematica
In[1]:= TrigFactor[Sin[a] Cos[b] + Cos[a] Sin[b]]
Out[1]= Sin[a + b]
```

```mathematica
In[1]:= TrigFactor[Cos[a] Cos[b] - Sin[a] Sin[b]]
Out[1]= Cos[a + b]
```

A difference of squares of sine and cosine folds into a double angle:

```mathematica
In[1]:= TrigFactor[Sin[x]^2 - Cos[x]^2]
Out[1]= -Cos[2 x]
```

A full trigonometric perfect square is factored and phase-shifted into a single
squared sine:

```mathematica
In[1]:= TrigFactor[Sin[x]^2 + 2 Sin[x] Cos[x] + Cos[x]^2]
Out[1]= 2 Sin[1/4 Pi + x]^2
```

It handles hyperbolic functions too, recognising the `cosh` double-angle identity:

```mathematica
In[1]:= TrigFactor[Sinh[x]^2 + Cosh[x]^2]
Out[1]= Cosh[2 x]
```
