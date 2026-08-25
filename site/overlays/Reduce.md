---
status: Partial
references:
  - "Collins & Hong, \"Partial Cylindrical Algebraic Decomposition for Quantifier Elimination\", J. Symbolic Computation 12 (1991)."
  - "McCallum, \"An Improved Projection Operation for Cylindrical Algebraic Decomposition\", in Quantifier Elimination and Cylindrical Algebraic Decomposition (1998)."
  - "Basu, Pollack & Roy, \"Algorithms in Real Algebraic Geometry\" (2nd ed., 2006), Ch. 5 & 11 (sign determination, CAD)."
  - "Schrijver, \"Theory of Linear and Integer Programming\" (1986), §12.2 (Fourier-Motzkin elimination)."
---
### Worked examples

Equations over the **complexes** (the default): `Reduce` keeps *every* case,
including the degenerate branches `Solve` drops.

```mathematica
In[1]:= Reduce[a x == b, x]
Out[1]= a != 0 && x == b/a || a == 0 && b == 0

In[2]:= Solve[a x == b, x]
Out[2]= {{x -> b/a}}
```

```mathematica
In[1]:= Reduce[x^2 == 4, x]
Out[1]= x == -2 || x == 2

In[2]:= Reduce[x^2 == -1, x]
Out[2]= x == -I || x == I

In[3]:= Reduce[a x^2 + b x + c == 0, x]
Out[3]= a != 0 && x == (1/2 (-b + Sqrt[b^2 - 4 a c]))/a || a != 0 && x == (1/2 (-b - Sqrt[b^2 - 4 a c]))/a || a == 0 && b != 0 && x == -c/b || a == 0 && b == 0 && c == 0
```

Inequalities over the **reals**. A statement with an ordering inequality and no
explicit domain defaults to `Reals`, so `Reals` may be omitted:

```mathematica
In[1]:= Reduce[x^2 > 1, x, Reals]
Out[1]= x < -1 || x > 1

In[2]:= Reduce[x^2 < 2, x, Reals]
Out[2]= -Sqrt[2] < x < Sqrt[2]

In[3]:= Reduce[(x - 1) (x - 2) (x - 3) > 0, x, Reals]
Out[3]= 1 < x < 2 || x > 3
```

A rational-function inequality treats the denominator roots as excluded poles,
and a chained inequality is solved directly:

```mathematica
In[1]:= Reduce[1/x < 1, x, Reals]
Out[1]= x < 0 || x > 1

In[2]:= Reduce[-5 < 3 x + 7 <= 22, x]
Out[2]= -4 < x <= 5

In[3]:= Reduce[Abs[x] < 1, x, Reals]
Out[3]= -1 < x < 1
```

A statement that decides returns a plain `True` or `False`:

```mathematica
In[1]:= Reduce[x^2 + 1 > 0, x, Reals]
Out[1]= True

In[2]:= Reduce[x < 0 && x > 1, x, Reals]
Out[2]= False
```

**Several variables.** Linear systems go through Fourier-Motzkin elimination;
genuinely nonlinear ones through a Cylindrical Algebraic Decomposition:

```mathematica
In[1]:= Reduce[x + y < 1 && x > 0 && y > 0, {x, y}, Reals]
Out[1]= 0 < x < 1 && 0 < y < 1 - x

In[2]:= Reduce[x y > 0, {x, y}, Reals]
Out[2]= x < 0 && y < 0 || x > 0 && y > 0
```

**Over the integers.** The `Integers`/`Rationals` domain reuses the `Solve`
Diophantine engine and reports its answer in logical form, introducing an
integer parameter `C[k]` for an infinite family:

```mathematica
In[1]:= Reduce[2 x + 3 y == 1, {x, y}, Integers]
Out[1]= Element[C[1], Integers] && x == -1 + 3 C[1] && y == 1 - 2 C[1]

In[2]:= Reduce[x^2 < 10 && x > 0, x, Integers]
Out[2]= x == 1 || x == 2 || x == 3
```

### Notes

`Reduce` and `Solve` answer different questions. `Solve[eqns, vars]` returns the
*generic* solution as a list of replacement rules and quietly discards the
special cases — `Solve[a x == b, x]` gives `{{x -> b/a}}`, silently assuming
`a != 0`. `Reduce[expr, vars]` returns instead a quantifier-free `And`/`Or` tree
of equations and inequalities that describes the **entire** solution set, every
parametric and boundary case included, so it reports both `a != 0 && x == b/a`
and the degenerate branch `a == 0 && b == 0`. Use `Solve` when you want to
substitute a solution back with `/.`; use `Reduce` when you need the complete
condition under which a statement holds — especially with symbolic parameters,
or with inequalities, where there is no finite list of rules to return.

**Domain.** The default domain is `Complexes`. Because an ordering (`<`, `<=`,
`>`, `>=`) is undefined over the complex numbers, a statement that contains one
and is given no explicit domain is solved over the `Reals` — matching
Mathematica — so the `Reals` argument is optional whenever the statement
already contains an inequality. Equations (`==`) and `Unequal` (`!=`) keep the
`Complexes` default. Pass a third argument (`Complexes`, `Reals`, `Integers`, or
`Rationals`) to choose explicitly.

**Soundness over completeness.** `Reduce` never returns a wrong or incomplete
formula: when a sign or ordering cannot be decided exactly, when it meets a
construct it does not yet handle, or when it reaches an engine not yet wired
(nonlinear equations over the complexes, and general quantifier elimination),
it leaves the input unevaluated rather than guess. An empty logical result
(`False`) is therefore a genuine proof that no solution exists, and a returned
formula is guaranteed to describe the whole solution set exactly.
