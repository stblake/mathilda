### Worked examples

```mathematica
In[1]:= FullForm[a + b]
Out[1]= Plus[a, b]
```

```mathematica
In[1]:= FullForm[1/2]
Out[1]= Rational[1, 2]
```

```mathematica
In[1]:= FullForm[x^2 + 1]
Out[1]= Plus[1, Power[x, 2]]
```

```mathematica
In[1]:= FullForm[a/b]
Out[1]= Times[a, Power[b, -1]]
```

```mathematica
In[1]:= FullForm[x_Integer /; x > 0]
Out[1]= Condition[Pattern[x, Blank[Integer]], Greater[x, 0]]
```

### Notes

`FullForm` reveals the raw internal tree, with every head written before its arguments and no special-cased syntax. It is the quickest way to see how surface notation like `+`, `/`, and `^` maps onto the underlying `Plus`/`Rational`/`Power` heads. Division is `Times` with a `Power[_, -1]` factor, so `a/b` is `Times[a, Power[b, -1]]`. It also exposes the desugaring of pattern syntax — `x_Integer /; x > 0` is really `Condition[Pattern[x, Blank[Integer]], Greater[x, 0]]` — which is invaluable when debugging why a rule does or does not match.
