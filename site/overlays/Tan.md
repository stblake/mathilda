### Worked examples

```mathematica
In[1]:= Tan[Pi/4]
Out[1]= 1

In[2]:= N[Tan[1]]
Out[2]= 1.55741

In[3]:= Tan[Pi/2]
Out[3]= ComplexInfinity
```

### Notes

`Tan[z]` is equivalent to `Sin[z]/Cos[z]`. Singularities at `z = Pi/2 + k Pi` yield `ComplexInfinity`. `Tan` is Listable.
