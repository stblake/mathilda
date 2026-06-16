### Worked examples

```mathematica
In[1]:= Level[a + b c, {1}]
Out[1]= {a, b c}
```

Level `{2}` reaches one step deeper, into the factors of `b c`:

```mathematica
In[1]:= Level[a + b c, {2}]
Out[1]= {b, c}
```

`{-1}` extracts the atomic leaves; `{-2}` the parts of depth two:

```mathematica
In[1]:= Level[f[g[h[x]]], {-1}]
Out[1]= {x}

In[2]:= Level[(1 + x)^2 + y, {-2}]
Out[2]= {1 + x}
```

A full depth-first walk traverses subexpressions in lexicographic index order:

```mathematica
In[1]:= Level[{{a, b}, {c, {d, e}}}, Infinity]
Out[1]= {a, b, {a, b}, c, d, e, {d, e}, {c, {d, e}}}
```

With `Heads -> True`, the head of each expression and its parts are included:

```mathematica
In[1]:= Level[f[g[x], y], 2, Heads -> True]
Out[1]= {f, g, x, g[x], y}
```

### Notes

`Level[expr, levelspec]` returns the subexpressions of `expr` on the specified
levels. A positive level `n` selects parts reachable by `n` indices; a negative
level `-n` selects parts of depth `n`; `{-1}` gives all atoms and level `0` is
the whole expression. `Level[expr, levelspec, f]` applies `f` to the sequence of
subexpressions, and `Heads -> True` includes expression heads in the traversal.
