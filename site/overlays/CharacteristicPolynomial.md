---
status: Stable
references:
  - "G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — the characteristic polynomial and the eigenvalue problem."
  - "D. K. Faddeev and V. N. Faddeeva, *Computational Methods of Linear Algebra*, W. H. Freeman, 1963 — the Faddeev–LeVerrier–Souriau recurrence."
---
### Worked examples

```mathematica
In[1]:= CharacteristicPolynomial[{{1, 2}, {3, 4}}, x]
Out[1]= -2 - 5 x + x^2
```

```mathematica
In[1]:= CharacteristicPolynomial[{{a, b}, {c, d}}, x]
Out[1]= -b c + a d - a x - d x + x^2
```

```mathematica
In[1]:= CharacteristicPolynomial[IdentityMatrix[3], x]
Out[1]= 1 - 3 x + 3 x^2 - x^3
```

```mathematica
In[1]:= CharacteristicPolynomial[{{1/3, 1/2, 3/5}, {1/2, 4/5, 1}, {3/5, 1, 9/7}}, x]
Out[1]= 1/10500 - 239/2100 x + 254/105 x^2 - x^3
```

```mathematica
In[1]:= CharacteristicPolynomial[{{{1, 2}, {5, 4}}, {{4, 3}, {6, 4}}}, x]
Out[1]= -6 + 7 x - 2 x^2
```

```mathematica
In[1]:= CharacteristicPolynomial[{{{1, 1, 1}, {1, 0, 1}, {0, 0, 1}}, {{0, 1, 1}, {0, 1, 1}, {1, 0, 0}}}, x]
Out[1]= -1 - x + x^2
```

### Notes

`CharacteristicPolynomial[m, x]` is `Det[m - x I]` and `CharacteristicPolynomial[{m, a}, x]` is `Det[m - x a]`, returned as an expanded polynomial in `x`. It is the polynomial whose roots are the (generalised) eigenvalues of `m`, so it pairs naturally with `Eigenvalues`.

The ordinary case is computed by the Faddeev–LeVerrier recurrence in `O(n^4)`, so the characteristic polynomial of a large numeric matrix is found in polynomial time rather than through an `O(n!)` symbolic determinant. For an odd-order matrix the polynomial is monic-negative (leading term `-x^n`), matching `Det[m - x I]`. In the generalised case a shared null space of `m` and `a` lowers the degree — the missing leading term corresponds to an infinite generalised eigenvalue, so `Eigenvalues[{m, a}]` returns `Infinity` for each degree drop. The second argument may also be a number or an expression, in which case the polynomial is evaluated at that value.
