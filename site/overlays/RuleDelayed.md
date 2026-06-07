### Worked examples

```mathematica
In[1]:= {1, 2, 3} /. n_ :> n^2
Out[1]= {1, 4, 9}

In[2]:= FullForm[a :> b]
Out[2]= RuleDelayed[a, b]
```

### Notes

`a :> b` is shorthand for `RuleDelayed[a, b]`. The right-hand side is held and evaluated separately for each match, after the pattern bindings are substituted in.
