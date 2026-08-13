# Variables

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Variables[poly]`**

gives the sorted list of independent variables that appear as bases of non-numeric subexpressions in poly.

<details>
<summary>Notes</summary>

Walks the expression tree and collects symbols and compound forms (e.g. Sin\[x\], a\[i\]) that occur outside numeric arithmetic; duplicates are removed via canonical order.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Variables[(x + y)^2 + 3 z^2 - y z + 7]
Out[1]= {x, y, z}

In[2]:= Variables[Sin[x] + Cos[x]]
Out[2]= {Cos[x], Sin[x]}

In[3]:= Variables[E^x]
Out[3]= {}
```

### Applications (5)

```mathematica
In[4]:= Variables[x^2 + y z]
Out[4]= {x, y, z}

In[5]:= Variables[a x^2 + b x + c]
Out[5]= {a, b, c, x}

In[6]:= Variables[Sin[x] + y]
Out[6]= {Sin[x], y}

In[7]:= Variables[x^2 + 3 x + 2]
Out[7]= {x}

In[8]:= Variables[x^(1/2) + y^(2/3)]
Out[8]= {x, y}
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

## References

**See also:** [Plus](../../arithmetic/Plus/), [Times](../../arithmetic/Times/), [Power](../../arithmetic/Power/), [List](../../other-advanced/List/), [Pi](../../mathematical-constants/Pi/), [E](../../mathematical-constants/E/), [I](../../mathematical-constants/I/)

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 3 (multivariate polynomial representation and variable sets).
- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_div2.c`](https://github.com/stblake/mathilda/blob/main/tests/test_div2.c)
- Tests: [`tests/test_integrate_goursat.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_goursat.c)
- Tests: [`tests/test_poly.c`](https://github.com/stblake/mathilda/blob/main/tests/test_poly.c)
- Tests: [`tests/test_simplify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_simplify.c)

## Notes & additional examples

### Notes

`Variables` collects the independent generators that appear as bases of
non-numeric subexpressions and returns them in canonical sorted order with
duplicates removed. Pure numeric coefficients are ignored, so
`x^2 + 3 x + 2` reports only `{x}`. Compound non-atomic forms are treated as
single generators rather than being broken open: `Sin[x] + y` yields
`{Sin[x], y}`, keeping `Sin[x]` whole. Every symbol that occurs outside
numeric arithmetic is included, which is why parameters like `a, b, c` appear
alongside the main variable `x`.
