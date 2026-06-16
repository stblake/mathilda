### Worked examples

```mathematica
In[1]:= Sec[Pi/3]
Out[1]= 2

In[2]:= Sec[0]
Out[2]= 1

In[3]:= N[Sec[1]]
Out[3]= 1.85082
```

Exact special values come out in closed form, including the golden-ratio-related
`Sec[Pi/5]`:

```mathematica
In[1]:= Sec[Pi/4]
Out[1]= Sqrt[2]

In[2]:= Sec[Pi/5]
Out[2]= -1 + Sqrt[5]
```

An imaginary argument folds onto the hyperbolic secant via `Sec[I z] = Sech[z]`:

```mathematica
In[1]:= Sec[I]
Out[1]= Sech[1]
```

The Maclaurin series of `Sec` exposes the secant (Euler) numbers `1, 5, 61, ...`
in its coefficients:

```mathematica
In[1]:= Series[Sec[x], {x, 0, 6}]
Out[1]= 1 + 1/2 x^2 + 5/24 x^4 + 61/720 x^6 + O[x]^7
```

High-precision evaluation is available through `N`:

```mathematica
In[1]:= N[Sec[1], 40]
Out[1]= 1.8508157176809256179117532413986501934704
```

### Notes

`Sec[z]` is `1/Cos[z]`. Singularities at `z = Pi/2 + k Pi` yield `ComplexInfinity`. `Sec` is Listable.
