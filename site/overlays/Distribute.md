### Worked examples

```mathematica
In[1]:= Distribute[(a + b)(c + d)]
Out[1]= a c + b c + a d + b d
```

```mathematica
In[1]:= Distribute[f[a + b, c + d]]
Out[1]= f[a, c] + f[b, c] + f[a, d] + f[b, d]
```

```mathematica
In[1]:= Distribute[And[a, Or[b, c]], Or, And]
Out[1]= a && b || a && c
```

```mathematica
In[1]:= Distribute[f[a + b + c], Plus, f, Plus, g]
Out[1]= g[a] + g[b] + g[c]
```

### Notes

`Distribute[f[x1, x2, ...]]` forms the sum of `f` applied to every term in the
Cartesian product of the summands found in the `xi`. On `Times` it reproduces
ordinary expansion (`(a+b)(c+d)`), but the mechanism is general: `f[a+b, c+d]`
distributes an arbitrary head over its `Plus` arguments. The two extra-argument
form `Distribute[expr, g]` distributes over the head `g` instead of `Plus`; the
three-argument form `Distribute[expr, g, f]` only acts when the head of `expr` is
`f`. The third example distributes `And` over `Or` — the logical distributive law,
turning a clause into disjunctive form. The five-argument form
`Distribute[expr, g, f, gp, fp]` substitutes `gp` and `fp` for `g` and `f` in the
result, here rewriting the distributed terms under a new function head `g`.
