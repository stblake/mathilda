### Worked examples

```mathematica
In[1]:= x /. x -> 2
Out[1]= 2

In[2]:= FullForm[a -> b]
Out[2]= Rule[a, b]

In[3]:= {1, 2, 3} /. n_Integer -> 0
Out[3]= {0, 0, 0}
```

The right-hand side of `->` is evaluated immediately, when the rule object is
built. Contrast `x -> r` (which captures the current value of `r`) with the
delayed `x :> r` (which re-reads `r` each time it fires):

```mathematica
In[1]:= r = 5; {x :> r, x -> r}
Out[1]= {x :> r, x -> 5}
```

Pairing a `Rule` with a symbolic computation evaluates a closed-form result at a
point. Here `D[...]` is taken symbolically, then `x -> 5` substitutes the value:

```mathematica
In[1]:= D[x^3 + 2 x, x] /. x -> 5
Out[1]= 77
```

### Notes

`a -> b` is shorthand for `Rule[a, b]`. The right-hand side of a `Rule` is evaluated immediately when the rule object is built, unlike `RuleDelayed` (`:>`).
