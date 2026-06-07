### Worked examples

```mathematica
In[1]:= FullForm[a + b]
Out[1]= Plus[a, b]

In[2]:= FullForm[1/2]
Out[2]= Rational[1, 2]

In[3]:= FullForm[x^2 + 1]
Out[3]= Plus[1, Power[x, 2]]
```

### Notes

`FullForm` reveals the raw internal tree, with every head written before its arguments and no special-cased syntax. It is the quickest way to see how surface notation like `+`, `/`, and `^` maps onto the underlying `Plus`/`Rational`/`Power` heads.
