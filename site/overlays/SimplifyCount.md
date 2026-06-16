### Worked examples

```mathematica
In[1]:= SimplifyCount[1]
Out[1]= 1

In[2]:= SimplifyCount[x + 1]
Out[2]= 3
```

```mathematica
In[1]:= SimplifyCount[3 x]
Out[1]= 3

In[2]:= SimplifyCount[12345]
Out[2]= 5
```

```mathematica
In[1]:= SimplifyCount[a b c + d]
Out[1]= 6

In[2]:= SimplifyCount[Sin[x]^2 + Cos[x]^2]
Out[1]= 9
```

### Notes

`SimplifyCount` is the default complexity measure that `Simplify` minimises when no `ComplexityFunction` option is supplied. It counts subexpressions: each leaf and each operator head contributes, integers contribute their decimal digit count plus a constant for the sign (so `12345` costs `5`), and reals contribute `2`. Because `Simplify` keeps whichever candidate transform yields the smallest count, this is the same yardstick that decides, for example, that `1 + x` beats `(x^2 - 1)/(x - 1)`.
