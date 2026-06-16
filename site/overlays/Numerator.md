---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), on rational normal forms."
---
### Worked examples

```mathematica
In[1]:= Numerator[6/8]
Out[1]= 3
```

```mathematica
In[1]:= Numerator[(x+1)/(x-1)]
Out[1]= 1 + x
```

```mathematica
In[1]:= Numerator[x^(-2)]
Out[1]= 1
```

```mathematica
In[1]:= Numerator[a/b + c/d]
Out[1]= a/b + c/d
```

Combine the sum into a single fraction first, then `Numerator` returns the
genuine combined top:

```mathematica
In[1]:= Numerator[Together[a/b + c/d]]
Out[1]= b c + a d
```

Several factors are sorted into numerator and denominator by the sign of their
exponents — only `z^(-1)` is pushed below the bar:

```mathematica
In[1]:= Numerator[2 x/(3 y) * z^(-1)]
Out[1]= 2 x
```

### Notes

Numerator extracts the numerator of the structural rational form of its argument.
A rational constant is first reduced to lowest terms, so `Numerator[6/8] = 3`. For
symbolic quotients it returns the literal top of the `expr/expr` form, giving
`1 + x` for `(x+1)/(x-1)`. Factors with negative exponents are treated as
denominators, so `Numerator[x^(-2)] = 1`. Note that Mathilda does not auto-combine
a sum into a single fraction first: `Numerator[a/b + c/d]` returns the unevaluated
sum, since the expression is a `Plus`, not a single quotient. Apply `Together`
first if you want the combined numerator.
