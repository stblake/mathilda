# CoefficientList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
CoefficientList[poly, var] gives a list of coefficients of powers of var in poly, starting with power 0.
CoefficientList[poly, {var1, var2, ...}] gives an array of coefficients of the variables.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= CoefficientList[1 + 6 x - x^4, x]
Out[1]= {1, 6, 0, 0, -1}

In[2]:= CoefficientList[(1 + x)^10, x]
Out[2]= {1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1}

In[3]:= CoefficientList[1 + a x^2 + b x y + c y^2, {x, y}]
Out[3]= {{1, 0, c}, {0, b, 0}, {a, 0, 0}}
```

## Implementation notes

- `Protected`.
- Gives an array of coefficients starting with power 0.
- Returns a full rectangular array for multiple variables. Combinations of powers that do not appear in `poly` give zeros in the array.
- Automatically expands the polynomial internally.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 3 (dense polynomial coefficient vectors).
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

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

### Notes

`CoefficientList[poly, var]` returns the dense coefficient vector ordered by
ascending power, starting at the constant term (power 0); missing powers are
filled with explicit zeros, as the `{1, 0, 0, 1}` for `1 + x^3` shows. The
list length is one more than the degree. With a list of variables, the result
is a nested array whose `[[i, j]]` entry is the coefficient of
`x^(i-1) y^(j-1)`, so `x^2 + x y + y^2` produces a 3x3 matrix with ones on the
anti-diagonal. Coefficients may be symbolic, preserving parameters such as
`a, b, c`.
