---
status: Stable
references:
---
### Worked examples

```mathematica
In[1]:= Nest[f, x, 3]
Out[1]= f[f[f[x]]]
```

```mathematica
In[1]:= Nest[#^2 &, 2, 3]
Out[1]= 256
```

```mathematica
In[1]:= Nest[1/(1 + #) &, x, 2]
Out[1]= 1/(1 + 1/(1 + x))
```

### Notes

`Nest[f, expr, n]` applies `f` to `expr` exactly `n` times and returns only the
final result. `n` must be a non-negative integer; `Nest[f, expr, 0]` returns
`expr` unchanged. Each intermediate application is evaluated before the next, so
numeric iterations like `#^2 &` collapse to a single number (`2 -> 4 -> 16 ->
256`). Use `NestList` instead when the intermediate values are also wanted.
