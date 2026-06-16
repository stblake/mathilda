### Worked examples

```mathematica
In[1]:= SolveAlways[a x^2 + b x + c == 0, x]
Out[1]= {{c -> 0, b -> 0, a -> 0}}
```

For the quadratic to vanish for *every* `x`, all three coefficients must be zero.

```mathematica
In[1]:= SolveAlways[(a + b) x + (a - b - 2) == 0, x]
Out[1]= {{a -> 1, b -> -1}}
```

```mathematica
In[1]:= SolveAlways[p x^2 + q x + r == (x - 1)(x - 2), x]
Out[1]= {{r -> 2, q -> -3, p -> 1}}
```

Matching `p x^2 + q x + r` to `(x-1)(x-2) = x^2 - 3x + 2` recovers the coefficients by equating powers of `x`.

### Notes

`SolveAlways[eqns, vars]` finds the parameters (symbols in `eqns` but not in `vars`) for which the equations hold identically in `vars`. Each `lhs - rhs` is expanded as a polynomial in `vars` via `CoefficientList`; every coefficient is required to vanish, and the resulting system is handed to `Solve` for the parameters. It returns `{}` when there are no parameters to solve for.
