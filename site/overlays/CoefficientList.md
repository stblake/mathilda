---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), Ch. 3 (dense polynomial coefficient vectors)."
---
### Worked examples

```mathematica
In[1]:= CoefficientList[x^2 + 3 x + 2, x]
Out[1]= {2, 3, 1}
```

```mathematica
In[1]:= CoefficientList[1 + x^3, x]
Out[1]= {1, 0, 0, 1}
```

```mathematica
In[1]:= CoefficientList[a + b x + c x^2, x]
Out[1]= {a, b, c}
```

```mathematica
In[1]:= CoefficientList[x^2 + x y + y^2, {x, y}]
Out[1]= {{0, 0, 1}, {0, 1, 0}, {1, 0, 0}}
```

```mathematica
In[1]:= CoefficientList[(1 + x)^6, x]
Out[1]= {1, 6, 15, 20, 15, 6, 1}
```

```mathematica
In[1]:= CoefficientList[Expand[(1 + x + x^2)^4], x]
Out[1]= {1, 4, 10, 16, 19, 16, 10, 4, 1}
```

### Notes

`CoefficientList[poly, var]` returns the dense coefficient vector ordered by
ascending power, starting at the constant term (power 0); missing powers are
filled with explicit zeros, as the `{1, 0, 0, 1}` for `1 + x^3` shows. The
list length is one more than the degree. With a list of variables, the result
is a nested array whose `[[i, j]]` entry is the coefficient of
`x^(i-1) y^(j-1)`, so `x^2 + x y + y^2` produces a 3x3 matrix with ones on the
anti-diagonal. Coefficients may be symbolic, preserving parameters such as
`a, b, c`. Applied to `(1 + x)^6` it reads off a row of Pascal's triangle
(`{1, 6, 15, 20, 15, 6, 1}`), and the (pre-expanded) trinomial `(1 + x + x^2)^4`
yields the symmetric central-trinomial coefficient row
`{1, 4, 10, 16, 19, 16, 10, 4, 1}` — a quick way to read off such combinatorial
sequences.
