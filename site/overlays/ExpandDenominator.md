### Worked examples

```mathematica
In[1]:= ExpandDenominator[(a + b)^2 / (c + d)^2]
Out[1]= (a + b)^2/(c^2 + 2 c d + d^2)
```

```mathematica
In[1]:= ExpandDenominator[1 / ((x + 1)(x + 2)(x + 3))]
Out[1]= 1/(6 + 11 x + 6 x^2 + x^3)
```

```mathematica
In[1]:= ExpandDenominator[((x + 1)(x + 2)) / (y (y + 1))]
Out[1]= ((1 + x) (2 + x))/(y + y^2)
```

### Notes

`ExpandDenominator` multiplies out only the denominator (the negative integer
powers), leaving the numerator factored — the mirror image of
`ExpandNumerator`. The product of three linear factors above expands to the
cubic `6 + 11 x + 6 x^2 + x^3` in the denominator. It threads over lists,
equations, inequalities, and logic functions.
