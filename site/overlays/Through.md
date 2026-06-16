### Worked examples

```mathematica
In[1]:= Through[(f + g)[x]]
Out[1]= f[x] + g[x]
```

```mathematica
In[1]:= Through[{Sin, Cos, Tan}[Pi/4]]
Out[1]= {1/Sqrt[2], 1/Sqrt[2], 1}
```

```mathematica
In[1]:= f[x_] := x^2; g[x_] := x + 1; Through[(f + g)[3]]
Out[1]= 13
```

### Notes

`Through[p[f1, f2, ...][x]]` distributes the trailing arguments across each inner
operator, producing `p[f1[x], f2[x], ...]`. With `p = Plus` it turns a sum of
functions into a function of a point — so once `f` and `g` are defined,
`Through[(f + g)[3]]` evaluates `f[3] + g[3]`. Applying it to a list of heads
(`{Sin, Cos, Tan}`) evaluates all of them at the same argument at once.
`Through[expr, h]` distributes only when the outer head equals `h`.
