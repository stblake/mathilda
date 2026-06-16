### Worked examples

```mathematica
In[1]:= Collect[a x^2 + b x^2 + c x + x, x]
Out[1]= (1 + c) x + (a + b) x^2
```

```mathematica
In[1]:= Collect[(1 + x + y)^3, x]
Out[1]= 1 + x^3 + 3 y + 3 y^2 + y^3 + x (3 + 6 y + 3 y^2) + x^2 (3 + 3 y)
```

```mathematica
In[1]:= Collect[Expand[(1 + x)^4 (1 + y)^2], x, Factor]
Out[1]= (1 + y)^2 + 4 x (1 + y)^2 + 6 x^2 (1 + y)^2 + 4 x^3 (1 + y)^2 + x^4 (1 + y)^2
```

```mathematica
In[1]:= Collect[(x + y + z)^2, {x, y}]
Out[1]= x^2 + y^2 + 2 y z + z^2 + x (2 y + 2 z)
```

### Notes

`Collect[expr, x]` expands `expr` and gathers terms by power of `x`, leaving each coefficient free of `x`. A third argument applies a function (e.g. `Factor` or `Simplify`) to each coefficient before reassembly; a list of variables groups by each in turn.
