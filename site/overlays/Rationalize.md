### Worked examples

```mathematica
In[1]:= Rationalize[0.5]
Out[1]= 1/2

In[2]:= Rationalize[N[Pi], 10^-4]
Out[2]= 333/106

In[3]:= Rationalize[1.2 + 6.7 x]
Out[3]= 6/5 + 67/10 x
```

### Notes

`Rationalize[x]` finds a nearby rational with small denominator (within `c/q^2`, `c = 10^-4`); the two-argument form `Rationalize[x, dx]` returns the smallest-denominator rational within `dx` of `x`. It threads over compound expressions.
