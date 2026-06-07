### Worked examples

```mathematica
In[1]:= $Pre = Hold
Out[1]= Hold

In[2]:= 1 + 1
Out[2]= Hold[1 + 1]
```

### Notes

`$Pre`, if set, is applied to every input expression after parsing and before
standard evaluation. Using a head with `HoldAll` (such as `Hold`) lets `$Pre`
intercept the unevaluated input; otherwise the argument is evaluated first and
the effect is indistinguishable from `$Post`. Unset by default.
