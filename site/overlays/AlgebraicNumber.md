### Worked examples

```mathematica
In[1]:= AlgebraicNumber[Sqrt[2], {1, 2}]  (* represents 1 + 2 Sqrt[2] in Q(Sqrt[2]) *)
```

```mathematica
In[1]:= AlgebraicNumber[Sqrt[2], {3, 0}]  (* only the constant term: a rational, reduced out *)
```

```mathematica
In[1]:= N[AlgebraicNumber[Sqrt[2], {1, 2}]]  (* treated as a numeric quantity *)
```

```mathematica
In[1]:= AlgebraicNumber[Sqrt[2], {1, 2}] + AlgebraicNumber[Sqrt[2], {3, 4}]  (* arithmetic in the same field *)
```

```mathematica
In[1]:= AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2], {1, 2}], x]  (* recover the defining polynomial *)
```

### Notes

`AlgebraicNumber[theta, {c0, c1, …, cn}]` denotes `c0 + c1 theta + … + cn theta^n`
in the field `Q(theta)`. The object is automatically reduced so that `theta` is
an algebraic integer and the coefficient list has length equal to the degree of
`theta`'s minimal polynomial; an object representing a rational number collapses
to explicit rational form.

Objects in the same field combine under `+`, `*` and integer powers. The head
carries `NHoldAll`, so `N` reaches the dedicated numeric branch (giving the
value to any precision) rather than numericalising the stored generator and
coefficients in place. Requires FLINT for the canonicalisation.
