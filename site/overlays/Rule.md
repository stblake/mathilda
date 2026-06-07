### Worked examples

```mathematica
In[1]:= x /. x -> 2
Out[1]= 2

In[2]:= FullForm[a -> b]
Out[2]= Rule[a, b]

In[3]:= {1, 2, 3} /. n_Integer -> 0
Out[3]= {0, 0, 0}
```

### Notes

`a -> b` is shorthand for `Rule[a, b]`. The right-hand side of a `Rule` is evaluated immediately when the rule object is built, unlike `RuleDelayed` (`:>`).
