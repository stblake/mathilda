### Worked examples

```mathematica
In[1]:= MatchQ[{}, {RepeatedNull[1]}]
Out[1]= True

In[2]:= MatchQ[{1, 1}, {RepeatedNull[1]}]
Out[2]= True
```

This is exactly what makes a *trailing optional argument* expressible as a single
pattern: `f[x, RepeatedNull[y]]` matches both `f[x]` and `f[x, y, y]`:

```mathematica
In[1]:= MatchQ[f[x], f[x, RepeatedNull[y]]]
Out[1]= True

In[2]:= MatchQ[f[x, y, y], f[x, RepeatedNull[y]]]
Out[2]= True
```

A bound applies on top of the "zero allowed" base case, so `Cases` keeps the
empty and short lists alike:

```mathematica
In[1]:= Cases[{{a}, {a, a}, {a, a, a}}, {RepeatedNull[a, 2]}]
Out[1]= {{a}, {a, a}}
```

### Notes

`RepeatedNull[p]` (postfix `p...`) is like `Repeated` but matches *zero* or more occurrences, so an empty sequence also succeeds (Out[1]). Use it when the repeated part is optional.
