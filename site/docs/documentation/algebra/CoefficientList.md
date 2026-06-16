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

**Algorithm.** `builtin_coefficientlist` (in `src/poly/poly.c`) returns the dense coefficient array of a polynomial. It `expr_expand`s the input, computes each variable's degree with `get_degree_poly`, then the recursive worker `coeff_list_rec` builds a nested `List` whose shape mirrors the variable order. At each level it pulls all coefficients `c_0..c_d` of the current variable — preferring the bulk extractor `get_all_coeffs_expanded` (single pass over the expanded form) and falling back to `get_coeff_expanded` per degree — and recurses on each coefficient for the next variable.

**Data structures.** A `int* max_degrees` array sizes each axis; coefficients are produced as `Expr*` and assembled into nested `List` nodes. The bulk path avoids the naive (degree+1)-passes-per-level cost.

- `Protected`.
- Gives an array of coefficients starting with power 0.
- Returns a full rectangular array for multiple variables. Combinations of powers that do not appear in `poly` give zeros in the array.
- Automatically expands the polynomial internally.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 3 (dense polynomial coefficient vectors).
- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
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
