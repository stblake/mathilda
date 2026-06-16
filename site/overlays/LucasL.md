### Worked examples

```mathematica
In[1]:= LucasL[10]
Out[1]= 123
```

```mathematica
In[1]:= LucasL[5, x]
Out[1]= 5 x + 5 x^3 + x^5
```

```mathematica
In[1]:= LucasL[100]
Out[1]= 792070839848372253127
```

```mathematica
In[1]:= LucasL[-7]
Out[1]= -29
```

```mathematica
In[1]:= N[LucasL[20]/LucasL[19], 20]
Out[1]= 1.6180339887498948482
```

### Notes

`LucasL[n]` is the `n`th Lucas number `L_n` (`L_0 = 2`, `L_1 = 1`, `L_k = L_{k-1} + L_{k-2}`); `LucasL[n, x]` is the Lucas polynomial `L_n(x)`. Integer orders use GMP fast doubling for arbitrary size (so `LucasL[100]` is exact), and negative orders follow `L_{-n} = (-1)^n L_n`. Consecutive ratios `L_{n+1}/L_n` converge to `GoldenRatio` (last example). Inexact or complex orders use the closed form `phi^n + Cos[Pi n] phi^-n`. Listable; symbolic orders stay unevaluated.
