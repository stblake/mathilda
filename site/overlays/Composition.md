### Worked examples

```mathematica
In[1]:= Composition[f, g, h][x]
Out[1]= f[g[h[x]]]
```

```mathematica
In[1]:= Composition[Sqrt, Abs][-16]
Out[1]= 4
```

```mathematica
In[1]:= Composition[f, InverseFunction[f]]
Out[1]= Identity
```

### Notes

`Composition[f1, f2, ...]` represents a function that acts innermost-first when applied. It has attributes `Flat` and `OneIdentity`, and can be entered as `f1 @* f2 @* ...`. Compositions containing `Identity` or `InverseFunction[f]` automatically simplify.
