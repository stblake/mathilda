### Worked examples

```mathematica
In[1]:= NSolve[x^2 - 2 == 0, x]  (* a univariate polynomial, solved through NRoots *)
Out[1]= {{x -> -1.41421}, {x -> 1.41421}}
```

```mathematica
In[1]:= NSolve[x^3 - 1 == 0, x, Reals]  (* the Reals domain keeps only the real root *)
Out[1]= {{x -> 1.0}}
```

```mathematica
In[1]:= NSolve[{x + y == 3, x - y == 1}, {x, y}]  (* a linear system *)
Out[1]= {{x -> 2.0, y -> 1.0}}
```

```mathematica
In[1]:= NSolve[{x^2 + y^2 == 1, y == x}, {x, y}]  (* a zero-dimensional nonlinear system *)
Out[1]= {{x -> 0.707107, y -> 0.707107}, {x -> -0.707107, y -> -0.707107}}
```

### Notes

`NSolve[expr, vars]` returns numerical solutions as a list of replacement-rule
lists; `{}` means no solutions and `{{}}` the universal solution. `vars` may be a
single variable or a list, and `NSolve[{e1, e2, ...}, vars]` is the conjunction
`e1 && e2 && ...`. `NSolve[expr, vars, Reals]` restricts to real solutions; the
default domain is the complexes.

A univariate polynomial equation is handed to `NRoots` (roots repeated by
multiplicity), so `NRoots`'s `PrecisionGoal`/`AccuracyGoal` contract governs the
answer. A square, zero-dimensional polynomial **system** uses a Groebner-basis
multiplication-matrix (Moeller-Stetter) eigenvalue method: a Groebner basis gives
the quotient ring, and the eigenvalues/eigenvectors of the multiplication maps
yield each coordinate, every candidate verified against the original residuals.
`Method -> "Symbolic"` instead uses lexicographic elimination (solve the
univariate generator, back-substitute, recurse). Other equations fall back to
symbolic `Solve` then numericalisation, with a univariate `FindRoot` grid as a
last resort.

A working precision may be given as a trailing positional argument or via
`WorkingPrecision`; integer, real, and complex coefficients are handled at
machine and arbitrary precision. `MaxRoots`, `VerifySolutions`, and
`RandomSeeding` (the seed for the generic linear form) are also accepted.
