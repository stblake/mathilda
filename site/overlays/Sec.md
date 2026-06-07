### Worked examples

```mathematica
In[1]:= Sec[Pi/3]
Out[1]= 2

In[2]:= Sec[0]
Out[2]= 1

In[3]:= N[Sec[1]]
Out[3]= 1.85082
```

### Notes

`Sec[z]` is `1/Cos[z]`. Singularities at `z = Pi/2 + k Pi` yield `ComplexInfinity`. `Sec` is Listable.
