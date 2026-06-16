### Worked examples

```mathematica
In[1]:= Cos[Pi/3]
Out[1]= 1/2

In[2]:= Cos[{0, Pi/3, Pi}]
Out[2]= {1, 1/2, -1}

In[3]:= Cos[ArcCos[x]]
Out[3]= x
```

```mathematica
In[1]:= Cos[Pi/12]
Out[1]= 1/4 (Sqrt[2] + Sqrt[6])

In[2]:= Cos[Pi/5]
Out[2]= 1/4 (1 + Sqrt[5])
```

```mathematica
In[1]:= TrigExpand[Cos[a + b]]
Out[1]= Cos[a] Cos[b] - Sin[a] Sin[b]
```

```mathematica
In[1]:= Series[Cos[x], {x, 0, 8}]
Out[1]= 1 - 1/2 x^2 + 1/24 x^4 - 1/720 x^6 + 1/40320 x^8 + O[x]^9
```

```mathematica
In[1]:= N[Cos[1], 40]
Out[1]= 0.54030230586813971740093660744297660373228
```

### Notes

The argument is in radians; rational multiples of `Pi` reduce to exact values while numeric inputs route to libm / MPFR. `Cos` is Listable.
