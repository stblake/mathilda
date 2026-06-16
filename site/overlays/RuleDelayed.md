### Worked examples

```mathematica
In[1]:= {1, 2, 3} /. n_ :> n^2
Out[1]= {1, 4, 9}

In[2]:= FullForm[a :> b]
Out[2]= RuleDelayed[a, b]
```

Because the RHS is evaluated per match, a delayed rule can transform each matched
subexpression by a formula in its bindings. This one rewrites every monomial
`x^n` to its antiderivative `x^(n+1)/(n+1)` — termwise integration via one rule:

```mathematica
In[1]:= x^4 + 3 x^2 + 1 /. x^n_ :> x^(n+1)/(n+1)
Out[1]= 1 + x^3 + 1/5 x^5
```

Applying a function to each binding is the typical reason to prefer `:>` over
`->`, since the RHS must be re-evaluated for each distinct match:

```mathematica
In[1]:= {f[1], f[2], f[3]} /. f[n_] :> n!
Out[1]= {1, 2, 6}
```

### Notes

`a :> b` is shorthand for `RuleDelayed[a, b]`. The right-hand side is held and evaluated separately for each match, after the pattern bindings are substituted in.
