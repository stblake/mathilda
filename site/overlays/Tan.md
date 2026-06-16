### Worked examples

```mathematica
In[1]:= Tan[Pi/4]
Out[1]= 1

In[2]:= N[Tan[1]]
Out[2]= 1.55741

In[3]:= Tan[Pi/2]
Out[3]= ComplexInfinity
```

```mathematica
In[1]:= Tan[Pi/12]
Out[1]= 2 - Sqrt[3]
```

```mathematica
In[1]:= N[Tan[1], 40]
Out[1]= 1.5574077246549022305069748074583601730872
```

```mathematica
In[1]:= Simplify[Tan[ArcSin[x]]]
Out[1]= x/Sqrt[1 - x^2]
```

### Notes

`Tan[z]` is equivalent to `Sin[z]/Cos[z]`. Singularities at `z = Pi/2 + k Pi` yield `ComplexInfinity`. `Tan` is Listable.
