### Worked examples

```mathematica
In[1]:= ExpandNumerator[(a + b)^2 / (c + d)^2]
Out[1]= (a^2 + 2 a b + b^2)/(c + d)^2
```

```mathematica
In[1]:= ExpandNumerator[((x + 1)(x + 2)) / (y (y + 1))]
Out[1]= (2 + 3 x + x^2)/(y (1 + y))
```

```mathematica
In[1]:= ExpandNumerator[(1 + x)^3 / x^2 == (1 + y)^2 / y]
Out[1]= (1 + 3 x + 3 x^2 + x^3)/x^2 == (1 + 2 y + y^2)/y
```

### Notes

`ExpandNumerator` distributes products and positive integer powers in the
numerator only, leaving the denominator factored. Unlike `Expand`, it does not
split the fraction into separate terms. It threads over lists and relations, so
both sides of the equation above have their numerators expanded independently.
