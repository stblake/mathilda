### Worked examples

```mathematica
In[1]:= FactorTermsList[6 x^2 + 4 x]
Out[1]= {2, 2 x + 3 x^2}

In[2]:= FactorTermsList[2 x^2 + 4 x + 2]
Out[2]= {2, 1 + 2 x + x^2}
```

```mathematica
In[1]:= FactorTermsList[a b x^2 + a b c x, {x}]
Out[1]= {1, a b, c x + x^2}
```

```mathematica
In[1]:= FactorTermsList[12 x^3 y + 18 x^2 y, {x, y}]
Out[1]= {6, 1, y, 3 x^2 + 2 x^3}
```

### Notes

`FactorTermsList[poly]` returns `{overall numerical factor, polynomial with that factor removed}`. The numerical content (here 2) is pulled out, leaving the primitive part; it does not factor the remaining polynomial further. Given a variable list `{x1, ..., xn}`, it stratifies the polynomial: the first element is the numeric content, the second is the factor depending on none of the variables, and the later elements depend on progressively more of them. In the last example `12 x^3 y + 18 x^2 y` separates as `6 · 1 · y · (3 x^2 + 2 x^3)`, isolating the variable-free numeric content `6`, the `y`-only factor, and finally the `x`-dependent remainder.
