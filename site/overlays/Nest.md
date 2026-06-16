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

A symbolic continued-fraction approximant, built by nesting `1/(1+#)`:

```mathematica
In[1]:= Nest[1/(1 + #) &, x, 2]
Out[1]= 1/(1 + 1/(1 + x))
```

Newton's iteration for `Sqrt[2]` kept fully exact — four steps from `1` already
give a 12-figure rational approximant:

```mathematica
In[1]:= Nest[(# + 2/#)/2 &, 1, 4]
Out[1]= 665857/470832
```

Nesting a radical builds a finite "nested radical" tower:

```mathematica
In[1]:= Nest[Sqrt[1 + #] &, x, 2]
Out[1]= Sqrt[1 + Sqrt[1 + x]]
```

### Notes

`Nest[f, expr, n]` applies `f` to `expr` exactly `n` times and returns only the
final result. `n` must be a non-negative integer; `Nest[f, expr, 0]` returns
`expr` unchanged. Each intermediate application is evaluated before the next, so
numeric iterations like `#^2 &` collapse to a single number (`2 -> 4 -> 16 ->
256`). Use `NestList` instead when the intermediate values are also wanted.
