### Worked examples

```mathematica
In[1]:= PossibleZeroQ[(x - 1) (x + 1) - (x^2 - 1)]
Out[1]= True
```

```mathematica
In[1]:= PossibleZeroQ[x^2 + 1]
Out[1]= False
```

```mathematica
In[1]:= PossibleZeroQ[Sin[x]^2 + Cos[x]^2 - 1]
Out[1]= True
```

```mathematica
In[1]:= PossibleZeroQ[Sqrt[2] + Sqrt[3] - Sqrt[5 + 2 Sqrt[6]]]
Out[1]= True
```

```mathematica
In[1]:= PossibleZeroQ[Log[2] + Log[3] - Log[6]]
Out[1]= True
```

### Notes

`PossibleZeroQ[expr]` uses combined symbolic and numerical heuristics to decide
whether `expr` is identically zero. It sees through polynomial cancellation
(`(x-1)(x+1) - (x^2-1) = 0`), the Pythagorean identity
`Sin[x]^2 + Cos[x]^2 - 1`, the nested-radical denesting
`Sqrt[2] + Sqrt[3] = Sqrt[5 + 2 Sqrt[6]]`, and the logarithm law
`Log[2] + Log[3] = Log[6]`. A nonzero expression like `x^2 + 1` returns
`False`. Because deciding whether a closed-form expression is exactly zero is
undecidable in general, `PossibleZeroQ` is a fast but not infallible test:
`True` strongly suggests a zero and `False` rules one out for the cases it can
analyse.
