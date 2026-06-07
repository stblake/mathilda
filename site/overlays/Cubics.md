### Worked examples

```mathematica
In[1]:= Solve[x^3 + x + 1 == 0, x]
Out[1]= {{x -> Root[1 + #1 + #1^3 &, 1]}, {x -> Root[1 + #1 + #1^3 &, 2]}, {x -> Root[1 + #1 + #1^3 &, 3]}}

In[2]:= Solve[x^3 + x + 1 == 0, x, Cubics -> True]
Out[2]= {{x -> -1/3 ((1/2 (27 + 3 Sqrt[93]))^(1/3) - 3/(1/2 (27 + 3 Sqrt[93]))^(1/3))}, {x -> -1/3 ((1/2 (27 + 3 Sqrt[93]))^(1/3) (-1/2 + (1/2*I) Sqrt[3]) - 3 1/((1/2 (27 + 3 Sqrt[93]))^(1/3) (-1/2 + (1/2*I) Sqrt[3])))}, {x -> -1/3 ((1/2 (27 + 3 Sqrt[93]))^(1/3) (-1/2 + (1/2*I) Sqrt[3])^2 - 3 1/((1/2 (27 + 3 Sqrt[93]))^(1/3) (-1/2 + (1/2*I) Sqrt[3])^2))}}
```

### Notes

`Cubics` is a `Solve` option, not a function. With the default
`Cubics -> False`, irreducible cubics are returned as held `Root[]` objects;
setting `Cubics -> True` forces the explicit Cardano radical formulas, as the
contrast above shows.
