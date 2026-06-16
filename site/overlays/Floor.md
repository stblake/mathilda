### Worked examples

```mathematica
In[1]:= Floor[7/2]
Out[1]= 3
```

```mathematica
In[1]:= Floor[-2.3]
Out[1]= -3
```

```mathematica
In[1]:= Floor[17, 5]
Out[1]= 15
```

```mathematica
In[1]:= Floor[{2.7, -2.7, 5/2, 11/3}]
Out[1]= {2, -3, 2, 3}
```

```mathematica
In[1]:= Floor[N[Pi, 40] 10^20]
Out[1]= 314159265358979323846
```

### Notes

`Floor[x]` rounds toward `-Infinity`; the two-argument `Floor[x, a]` gives the
greatest multiple of `a` not exceeding `x`. Exact inputs return exact
integers, and `Floor` is `Listable` so it threads over a vector of mixed
reals and rationals. The last example extracts the first 21 significant
digits of `π` exactly by flooring a 40-digit MPFR value scaled by `10^20`.
