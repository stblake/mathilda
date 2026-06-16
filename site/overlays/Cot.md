### Worked examples

```mathematica
In[1]:= Cot[Pi/4]
Out[1]= 1
```

```mathematica
In[1]:= Cot[Pi/12]
Out[1]= 2 + Sqrt[3]
```

```mathematica
In[1]:= N[Cot[1], 40]
Out[1]= 0.64209261593433070300641998659426562023026
```

```mathematica
In[1]:= Series[Cot[x], {x, 0, 5}]
Out[1]= 1/x - 1/3 x - 1/45 x^3 - 2/945 x^5 + O[x]^6
```

```mathematica
In[1]:= Cot[I]
Out[1]= -I Coth[1]
```

### Notes

`Cot[z]` is equivalent to `Cos[z]/Sin[z]`. Singularities at `z = k Pi` yield `ComplexInfinity`. `Cot` is Listable. The Laurent series at the origin exposes the Bernoulli-number coefficients, and imaginary arguments fold onto the hyperbolic cotangent.
