### Worked examples

```mathematica
In[1]:= Csc[Pi/6]
Out[1]= 2

In[2]:= N[Csc[1]]
Out[2]= 1.1884

In[3]:= Csc[Pi]
Out[3]= ComplexInfinity
```

### Notes

`Csc[z]` is `1/Sin[z]`. Singularities at `z = k Pi` yield `ComplexInfinity`. `Csc` is Listable.
