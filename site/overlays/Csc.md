### Worked examples

```mathematica
In[1]:= Csc[Pi/6]
Out[1]= 2
```

```mathematica
In[1]:= Csc[Pi/12]
Out[1]= Sqrt[2] (1 + Sqrt[3])
```

```mathematica
In[1]:= N[Csc[1], 40]
Out[1]= 1.1883951057781212162615994523745510035279
```

```mathematica
In[1]:= Series[Csc[x], {x, 0, 5}]
Out[1]= 1/x + 1/6 x + 7/360 x^3 + 31/15120 x^5 + O[x]^6
```

```mathematica
In[1]:= Csc[I]
Out[1]= -I Csch[1]
```

### Notes

`Csc[z]` is `1/Sin[z]`. Singularities at `z = k Pi` yield `ComplexInfinity`. `Csc` is Listable. Exact special angles are returned in closed radical form, the Laurent expansion gives the cosecant's pole at the origin, and imaginary arguments map onto the hyperbolic cosecant.
