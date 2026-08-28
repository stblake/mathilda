---
status: Partial
references:
  - "Collins, \"Quantifier Elimination for Real Closed Fields by Cylindrical Algebraic Decomposition\", LNCS 33 (1975)."
  - "Collins & Hong, \"Partial Cylindrical Algebraic Decomposition for Quantifier Elimination\", J. Symbolic Computation 12 (1991)."
  - "McCallum, \"An Improved Projection Operation for Cylindrical Algebraic Decomposition\", in Caviness & Johnson (eds.), Quantifier Elimination and CAD (1998)."
  - "Basu, Pollack & Roy, \"Algorithms in Real Algebraic Geometry\" (2nd ed., 2006), Ch. 5 & 11."
---
### Worked examples

A single-variable inequality: the solution set is a union of intervals with exact endpoints.

```mathematica
In[1]:= CylindricalDecomposition[x^2 > 1, {x}]
Out[1]= x < -1 || x > 1

In[2]:= CylindricalDecomposition[x^3 - x > 0, {x}]
Out[2]= -1 < x < 0 || x > 1
```

With several variables the answer is *cylindrical*: each variable is bounded in terms of the ones before it. Here `y` is fenced by functions of `x`.

```mathematica
In[1]:= CylindricalDecomposition[x y > 1, {x, y}]
Out[1]= x < 0 && y < 1/x || x > 0 && y > 1/x
```

The closed unit disk, sliced along `x`:

```mathematica
In[1]:= CylindricalDecomposition[x^2 + y^2 <= 1, {x, y}]
Out[1]= -1 <= x <= 1 && -1/2 Sqrt[4 - 4 x^2] <= y <= 1/2 Sqrt[4 - 4 x^2]
```

When the statement is decided outright -- true everywhere, or nowhere -- the result collapses to a plain `True` or `False`:

```mathematica
In[1]:= CylindricalDecomposition[x^2 + 1 > 0, {x}]
Out[1]= True

In[2]:= CylindricalDecomposition[x^2 < 0, {x}]
Out[2]= False
```

Equations decompose to their real roots just as inequalities decompose to intervals:

```mathematica
In[1]:= CylindricalDecomposition[x^2 == 1, {x}]
Out[1]= x == -1 || x == 1
```

### Notes

`CylindricalDecomposition[expr, vars]` returns a *cylindrical algebraic decomposition*
of the real solution set of `expr` -- a quantifier-free `And`/`Or` formula in which the
variables are bounded one after another, each in terms of the earlier ones, so that the
region is described as a stack of cylinders. It is the geometric companion of
`Reduce` over the reals: the two return the same formula for the same real problem, and
`CylindricalDecomposition` exists to name the decomposition directly and to keep the
domain fixed at the `Reals` (an explicit `Reals` third argument is accepted but
redundant).

Because the domain is always the reals, the endpoints are real algebraic numbers and may
be irrational surds, as in the disk example above. A statement that decides returns
`True` or `False`; one whose decomposition cannot be computed exactly -- an undecidable
sign, or a positive-dimensional system with irrational fibres beyond the current engine
-- is returned unevaluated rather than guessed, so a formula that *is* returned describes
the whole real solution set exactly. All of `Reduce`'s options (`Modulus`, `Cubics`,
`Quartics`, `WorkingPrecision`, ...) may be given and are forwarded unchanged.
