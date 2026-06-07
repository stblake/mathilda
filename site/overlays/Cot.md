### Worked examples

```mathematica
In[1]:= Cot[Pi/4]
Out[1]= 1

In[2]:= N[Cot[1]]
Out[2]= 0.642093

In[3]:= Cot[Pi]
Out[3]= ComplexInfinity
```

### Notes

`Cot[z]` is equivalent to `Cos[z]/Sin[z]`. Singularities at `z = k Pi` yield `ComplexInfinity`. `Cot` is Listable.
