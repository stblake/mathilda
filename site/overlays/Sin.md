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

```mathematica
In[1]:= Sin[Pi/10]
Out[1]= 1/4 (-1 + Sqrt[5])

In[2]:= Sin[Pi/12]
Out[2]= 1/4 (Sqrt[6] - Sqrt[2])
```

```mathematica
In[1]:= Sin[I]
Out[1]= I Sinh[1]

In[2]:= TrigExpand[Sin[3 x]]
Out[2]= -Sin[x]^3 + 3 Cos[x]^2 Sin[x]

In[3]:= N[Sin[1], 40]
Out[3]= 0.84147098480789650665250232163029899962254
```

### Notes

The argument is in radians; rational multiples of `Pi` reduce to exact values while numeric inputs route to libm (Real) or MPFR. `Sin` is Listable, so it threads over lists element-wise.
