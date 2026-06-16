# Variables

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Variables[poly]
    gives the sorted list of independent variables that appear as bases
    of non-numeric subexpressions in poly.
Walks the expression tree and collects symbols and compound forms
(e.g. Sin[x], a[i]) that occur outside numeric arithmetic; duplicates
are removed via canonical order.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Variables[(x + y)^2 + 3 z^2 - y z + 7]
Out[1]= {x, y, z}

In[2]:= Variables[Sin[x] + Cos[x]]
Out[2]= {Cos[x], Sin[x]}

In[3]:= Variables[E^x]
Out[3]= {}
```

## Implementation notes

**Algorithm.** `builtin_variables` walks the expression with `collect_variables`, gathering the
distinct symbols that occur as polynomial generators (bare symbols and non-numeric bases,
excluding numeric atoms and known constants), then sorts the collected `Expr*` array with
`qsort` under `compare_expr_ptrs` (the canonical `expr_compare` order) and wraps the result in a
`List`. The output is the deduplicated, canonically-ordered list of variables on which the
input is treated as a polynomial/rational expression.

- `Protected`.
- Looks for variables only inside `Plus`, `Times`, and `Power` with rational exponents.
- Returns a sorted `List` of variables.
- Symbolic constants like `Pi`, `E`, and `I` are not treated as variables.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 3 (multivariate polynomial representation and variable sets).
- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Variables[x^2 + y z]
Out[1]= {x, y, z}
```

```mathematica
In[1]:= Variables[a x^2 + b x + c]
Out[1]= {a, b, c, x}
```

```mathematica
In[1]:= Variables[Sin[x] + y]
Out[1]= {Sin[x], y}
```

```mathematica
In[1]:= Variables[x^2 + 3 x + 2]
Out[1]= {x}
```

Fractional powers are treated as polynomial generators in their own right, so each radical base contributes its variable:

```mathematica
In[1]:= Variables[x^(1/2) + y^(2/3)]
Out[1]= {x, y}
```

### Notes

`Variables` collects the independent generators that appear as bases of
non-numeric subexpressions and returns them in canonical sorted order with
duplicates removed. Pure numeric coefficients are ignored, so
`x^2 + 3 x + 2` reports only `{x}`. Compound non-atomic forms are treated as
single generators rather than being broken open: `Sin[x] + y` yields
`{Sin[x], y}`, keeping `Sin[x]` whole. Every symbol that occurs outside
numeric arithmetic is included, which is why parameters like `a, b, c` appear
alongside the main variable `x`.
