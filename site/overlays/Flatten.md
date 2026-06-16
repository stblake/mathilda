### Worked examples

```mathematica
In[1]:= Flatten[{{1, 2}, {3, {4}}}]
Out[1]= {1, 2, 3, 4}
```

```mathematica
In[1]:= Flatten[{{1, {2, 3}}, {4, {5}}}, 1]
Out[1]= {1, {2, 3}, 4, {5}}
```

```mathematica
In[1]:= Flatten[Table[{i, j}, {i, 2}, {j, 2}], 1]
Out[1]= {{1, 1}, {1, 2}, {2, 1}, {2, 2}}
```

```mathematica
In[1]:= Flatten[f[a, f[b, f[c, d]]], 2, f]
Out[1]= f[a, b, c, d]
```

### Notes

Without a level argument, `Flatten` removes all nesting; with a level `n` it
flattens only the top `n` levels. The third example collapses a 2x2 nested
`Table` into a flat list of coordinate pairs (a common reshaping idiom). The
last shows the three-argument form `Flatten[expr, n, h]`, which flattens
nested calls of an arbitrary head `h` rather than `List`.
