### Worked examples

```mathematica
In[1]:= Solve[x^4 + x + 1 == 0, x]
Out[1]= {{x -> Root[1 + #1 + #1^4 &, 1]}, {x -> Root[1 + #1 + #1^4 &, 2]}, {x -> Root[1 + #1 + #1^4 &, 3]}, {x -> Root[1 + #1 + #1^4 &, 4]}}

In[2]:= Quartics
Out[2]= Quartics
```

A biquadratic quartic factors automatically into explicit radicals regardless of the option setting:

```mathematica
In[1]:= Solve[x^4 - 5 x^2 + 6 == 0, x]
Out[1]= {{x -> -Sqrt[2]}, {x -> Sqrt[2]}, {x -> -Sqrt[3]}, {x -> Sqrt[3]}}
```

A pure fourth power yields the four complex fourth roots in closed form:

```mathematica
In[1]:= Solve[x^4 - 2 == 0, x]
Out[1]= {{x -> -2^(1/4)}, {x -> 2^(1/4)}, {x -> -I 2^(1/4)}, {x -> I 2^(1/4)}}
```

A non-biquadratic but radical-solvable quartic gives nested radicals; here the roots are the four values `±Sqrt[2] ± Sqrt[3]` written as `Sqrt[(10 ± 4 Sqrt[6])/2]`:

```mathematica
In[1]:= Solve[x^4 - 10 x^2 + 1 == 0, x]
Out[1]= {{x -> -Sqrt[1/2 (10 - 4 Sqrt[6])]}, {x -> Sqrt[1/2 (10 - 4 Sqrt[6])]}, {x -> -Sqrt[1/2 (10 + 4 Sqrt[6])]}, {x -> Sqrt[1/2 (10 + 4 Sqrt[6])]}}
```

### Notes

`Quartics` is a `Solve` option (the quartic analogue of `Cubics`), not a
function; evaluating the bare symbol just returns itself. With the default
`Quartics -> False`, an irreducible quartic is returned as held `Root[]`
objects as above; `Quartics -> True` requests explicit radical formulas where
they apply (biquadratic and other special quartics still reduce automatically).
