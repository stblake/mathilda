### Worked examples

```mathematica
In[1]:= Tanh[0]
Out[1]= 0

In[2]:= N[Tanh[1]]
Out[2]= 0.761594

In[3]:= Tanh[ArcTanh[z]]
Out[3]= z
```

```mathematica
In[1]:= N[Tanh[1], 40]
Out[1]= 0.76159415595576488811945828260479359041279
```

```mathematica
In[1]:= D[Tanh[x], x]
Out[1]= Sech[x]^2
```

```mathematica
In[1]:= Series[Tanh[x], {x, 0, 7}]
Out[1]= x - 1/3 x^3 + 2/15 x^5 - 17/315 x^7 + O[x]^8
```

### Notes

`Tanh[z]` is the hyperbolic tangent, `Sinh[z]/Cosh[z]`. `Tanh` is Listable.
