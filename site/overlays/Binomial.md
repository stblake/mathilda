---
status: Stable
references:
  - "Knuth, \"The Art of Computer Programming, Vol. 2: Seminumerical Algorithms\", on binomial coefficients."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), on the generalized binomial and its polynomial form."
---
### Worked examples

```mathematica
In[1]:= Binomial[50, 25]
Out[1]= 126410606437752
```

```mathematica
In[1]:= Binomial[-1, 3]
Out[1]= -1
```

```mathematica
In[1]:= Binomial[1/2, 3]
Out[1]= 1/16
```

```mathematica
In[1]:= Binomial[n, 2]
Out[1]= 1/2 n (-1 + n)
```

### Notes

Integer arguments give exact coefficients via the falling-factorial product, with
`Binomial[50, 25]` returning the full bigint `126410606437752`. The generalized
definition extends to negative and rational upper arguments: `Binomial[-1, 3] =
-1` and `Binomial[1/2, 3] = 1/16`, using `binomial(x, k) = x(x-1)...(x-k+1)/k!`.
With a symbolic upper argument and a non-negative integer lower argument,
Binomial expands to a polynomial, so `Binomial[n, 2]` becomes `n(n-1)/2`.
