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

### Notes

`Solve` returns a list of solution rule-lists, one per solution; each inner
list assigns every requested variable. Complex roots are produced by default,
so `x^2 + 1 == 0` yields the conjugate pair `±I`, and irrational roots come
back in exact radical form (`±Sqrt[2]`). Linear systems are solved directly
and return a single rule-list. Cubic roots are reported using `(-1)^(1/3)`
style radicals; pass `Cubics -> False` / `Quartics -> False` to suppress
explicit radical forms when desired.
