---
status: Stable
references:
---
### Worked examples

```mathematica
In[1]:= Table[i^2, {i, 1, 5}]
Out[1]= {1, 4, 9, 16, 25}
```

```mathematica
In[1]:= Table[i + j, {i, 1, 2}, {j, 1, 3}]
Out[1]= {{2, 3, 4}, {3, 4, 5}}
```

```mathematica
In[1]:= Table[x, 4]
Out[1]= {x, x, x, x}
```

```mathematica
In[1]:= Table[i, {i, 0, 1, 1/2}]
Out[1]= {0, 1/2, 1}
```

```mathematica
In[1]:= Table[Fibonacci[n], {n, 1, 12}]
Out[1]= {1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144}
```

```mathematica
In[1]:= Table[Sum[1/k, {k, 1, n}], {n, 1, 5}]
Out[1]= {1, 3/2, 11/6, 25/12, 137/60}
```

```mathematica
In[1]:= Table[If[i == j, 1, 0], {i, 1, 3}, {j, 1, 3}]
Out[1]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}
```

### Notes

The single-argument iterator `{i, imax}` runs `i` from 1 to `imax`; `{i, imin,
imax}` and `{i, imin, imax, step}` give explicit bounds and a step. The step may
be an exact rational, so the values stay exact (`{0, 1/2, 1}`). Multiple iterator
specifications nest: the leftmost varies slowest, producing a list of lists.
`Table[expr, n]` with a plain count simply repeats `expr` `n` times. `Table` holds
its arguments, so the body is only evaluated as each iterator value is assigned —
the body may itself be a `Sum`, `D`, or any other computation (giving exact
harmonic numbers, identity matrices, and the like).
