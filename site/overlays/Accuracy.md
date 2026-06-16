### Worked examples

```mathematica
In[1]:= Accuracy[1.5]
Out[1]= 15.7785

In[2]:= Accuracy[5]
Out[2]= Infinity

In[3]:= Accuracy[N[Pi, 30]]
Out[3]= 29.6058
```

```mathematica
In[1]:= Accuracy[N[10^20, 30]]
Out[1]= 10.103

In[2]:= Accuracy[N[1/10^20, 30]]
Out[2]= 50.103
```

```mathematica
In[1]:= Accuracy[0.]
Out[1]= 323.607
```

### Notes

`Accuracy` gives the number of digits to the right of the decimal point, i.e. the digit count relative to the absolute (not relative) magnitude, satisfying `Accuracy[x] == Precision[x] - Log10[Abs[x]]`. Exact numbers have `Infinity` accuracy. For a fixed precision of 30 digits, the same `Pi`-style relative precision yields very different accuracies as magnitude changes: `N[10^20, 30]` has only `10.1` digits past the point, while `N[1/10^20, 30]` has `50.1`. Inexact zero is the special case — machine `0.` returns the finite value `323.607` rather than `Infinity`.
