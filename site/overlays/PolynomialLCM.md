### Worked examples

```mathematica
In[1]:= PolynomialLCM[x^2 - 1, x^2 + 2 x + 1]
Out[1]= (-1 + x) (1 + 2 x + x^2)
```

```mathematica
In[1]:= PolynomialLCM[x^6 - 1, x^4 - 1]
Out[1]= (-1 + x^6) (1 + x^2)
```

```mathematica
In[1]:= PolynomialLCM[x^2 - 2, x^2 - Sqrt[2], Extension -> Sqrt[2]]
Out[1]= 2 Sqrt[2] - 2 x^2 - Sqrt[2] x^2 + x^4
```

### Notes

`PolynomialLCM` returns the least common multiple of its polynomial arguments,
computed as `a b / gcd(a, b)`. For `x^2 - 1 = (x-1)(x+1)` and
`(x+1)^2` the shared factor `x+1` appears only once in the LCM, so the result
factors as `(x-1)(x+1)^2`. For `x^6-1` and `x^4-1` the common part is
`x^2-1`, leaving `(x^6-1)(x^2+1)`. The `Extension -> alpha` option computes
the LCM over `Q(alpha)`; with `alpha = Sqrt[2]` the inputs `x^2-2` and
`x^2-Sqrt[2]` are coprime there, so the LCM is their full product
`(x^2-2)(x^2-Sqrt[2])`.
