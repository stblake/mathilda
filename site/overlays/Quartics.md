### Worked examples

```mathematica
In[1]:= Solve[x^4 + x + 1 == 0, x]
Out[1]= {{x -> Root[1 + #1 + #1^4 &, 1]}, {x -> Root[1 + #1 + #1^4 &, 2]}, {x -> Root[1 + #1 + #1^4 &, 3]}, {x -> Root[1 + #1 + #1^4 &, 4]}}

In[2]:= Quartics
Out[2]= Quartics
```

### Notes

`Quartics` is a `Solve` option (the quartic analogue of `Cubics`), not a
function; evaluating the bare symbol just returns itself. With the default
`Quartics -> False`, an irreducible quartic is returned as held `Root[]`
objects as above; `Quartics -> True` requests explicit radical formulas where
they apply (biquadratic and other special quartics still reduce automatically).
