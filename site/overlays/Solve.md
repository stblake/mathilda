---
status: Stable
references:
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\" (3rd ed.), Ch. 14 (polynomial roots and resolution)."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), Ch. 9 (solving systems)."
---
### Worked examples

```mathematica
In[1]:= Solve[x^2 - 5 x + 6 == 0, x]
Out[1]= {{x -> 2}, {x -> 3}}
```

```mathematica
In[1]:= Solve[x^2 + 1 == 0, x]
Out[1]= {{x -> -I}, {x -> I}}
```

```mathematica
In[1]:= Solve[x^2 - 2 == 0, x]
Out[1]= {{x -> -Sqrt[2]}, {x -> Sqrt[2]}}
```

```mathematica
In[1]:= Solve[{x + y == 3, x - y == 1}, {x, y}]
Out[1]= {{x -> 2, y -> 1}}
```

```mathematica
In[1]:= Solve[x y == 1 && x + y == 3, {x, y}]
Out[1]= {{x -> 1/2 (3 - Sqrt[5]), y -> 1/2 (3 + Sqrt[5])}, {x -> 1/2 (3 + Sqrt[5]), y -> 1/2 (3 - Sqrt[5])}}

In[2]:= Solve[x y == 6 && x + y == 5, {x, y}, Integers]
Out[2]= {{x -> 3, y -> 2}, {x -> 2, y -> 3}}
```

```mathematica
In[1]:= Solve[a x^2 + b x + c == 0, x]
Out[1]= {{x -> (1/2 (-b + Sqrt[b^2 - 4 a c]))/a}, {x -> (1/2 (-b - Sqrt[b^2 - 4 a c]))/a}}
```

```mathematica
In[1]:= Solve[x^4 - 1 == 0, x]
Out[1]= {{x -> -1}, {x -> 1}, {x -> -I}, {x -> I}}

In[2]:= Solve[Sin[x] == 0, x]
Out[2]= {{x -> ConditionalExpression[Pi + 2 C[1] Pi, Element[C[1], Integers]]}, {x -> ConditionalExpression[2 C[1] Pi, Element[C[1], Integers]]}}
```

### Notes

`Solve` returns a list of solution rule-lists, one per solution; each inner
list assigns every requested variable. Complex roots are produced by default,
so `x^2 + 1 == 0` yields the conjugate pair `±I`, and irrational roots come
back in exact radical form (`±Sqrt[2]`). Linear systems are solved directly
and return a single rule-list. Cubic roots are reported using `(-1)^(1/3)`
style radicals; pass `Cubics -> False` / `Quartics -> False` to suppress
explicit radical forms when desired.

Nonlinear polynomial systems with finitely many solutions (a zero-dimensional
ideal) are solved via a lexicographic Gröbner basis and triangular
back-substitution, honouring the `Reals` / `Integers` domain. Systems with
infinitely many solutions (positive-dimensional ideals, e.g.
`Solve[x^2 - y^2 == 0, {x, y}]`) emit `Solve::nsdim` and are left unevaluated.
