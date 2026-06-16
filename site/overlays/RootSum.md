### Worked examples

```mathematica
In[1]:= RootSum[Function[t, t^2 + t + 1], Function[t, t^3]]
Out[1]= RootSum[Function[t, t^2 + t + 1], Function[t, t^3]]

In[2]:= RootSum[#^2 + 1 &, # &]
Out[2]= RootSum[#1^2 + 1 &, #1 &]
```

```mathematica
In[1]:= Integrate[1/(x^3 - x + 1), x]
Out[1]= RootSum[1 + #1^3 - #1 &, Log[-#1 + x]/(-1 + 3 #1^2) &]
```

```mathematica
In[1]:= D[RootSum[Function[t, t^3 + t + 1], Function[t, Log[x - t]/t]], x]
Out[1]= RootSum[Function[t, t^3 + t + 1], Function[t, 1/(t (-t + x))]]
```

### Notes

`RootSum[f, form]` denotes the formal sum of `form[a]` over the roots `a` of
`f[a] == 0` and is kept as a held symbolic object. It is produced by the
rational integrator's `NaiveLogPart` fallback when the logarithmic part of an
integral cannot be written with real elementary functions; differentiation
threads through the `form` function.
