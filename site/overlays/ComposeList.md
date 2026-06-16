### Worked examples

```mathematica
In[1]:= ComposeList[{f, g, h}, x]
Out[1]= {x, f[x], g[f[x]], h[g[f[x]]]}
```

```mathematica
In[1]:= ComposeList[{Sin, Cos, Tan}, x]
Out[1]= {x, Sin[x], Cos[Sin[x]], Tan[Cos[Sin[x]]]}
```

```mathematica
In[1]:= ComposeList[{1 + #1 &, #1^2 &, 2 #1 &}, a]
Out[1]= {a, 1 + a, (1 + a)^2, 2 (1 + a)^2}
```

### Notes

`ComposeList[{f1, f2, ...}, x]` builds `{x, f1[x], f2[f1[x]], ...}`, recording every intermediate of an innermost-first composition. The result has one more element than the list of functions, which makes it convenient for tracing iterated maps.
