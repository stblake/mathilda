# PolynomialQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PolynomialQ[expr, var]
    gives True if expr is a polynomial in var, False otherwise.
PolynomialQ[expr, {v1, v2, ...}]
    gives True if expr is a polynomial in all of the vi simultaneously.
Checks that expr expands to a sum of products of non-negative integer
powers of the vars with var-free coefficients.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PolynomialQ[x^3 - 2x/y + 3x z, x]
Out[1]= True

In[2]:= PolynomialQ[x^3 - 2x/y + 3x z, y]
Out[2]= False

In[3]:= PolynomialQ[x^2 + a x y^2 - b Sin[c], {x, y}]
Out[3]= True

In[4]:= PolynomialQ[f[a] + f[a]^2, f[a]]
Out[4]= True
```

## Implementation notes

**Algorithm.** `builtin_polynomialq` is a structural predicate. It normalises the second
argument into a variable set (a single symbol or a `List` of symbols) and calls
`is_polynomial`, which recurses over the expression tree: an expression is a polynomial in the
given variables iff every node is one of the variables, a sub-expression free of all the
variables (a degree-0 constant), a `Plus`/`Times` whose arguments are all polynomials, or a
`Power[base, k]` with a non-negative integer exponent `k` and polynomial base. Anything else
containing a variable in a non-polynomial position (e.g. `Sin[x]`, `1/x`, `x^(1/2)`) returns
`False`. Returns the symbol `True` or `False`.

- `Protected`.
- Variables can be symbols or compound expressions.
- Constants (expressions free of the specified variables) are polynomials of degree 0.
- `Power[base, exp]` is a polynomial if `exp` is a non-negative integer and `base` is a polynomial.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= PolynomialQ[x^3 + 2 x + 1, x]
Out[1]= True
```

```mathematica
In[1]:= PolynomialQ[Sin[x] + x, x]
Out[1]= False
```

```mathematica
In[1]:= PolynomialQ[x^2 y + x y^2 + 1, {x, y}]
Out[1]= True
```

```mathematica
In[1]:= PolynomialQ[x^2 + y/x, x]
Out[1]= False
```

### Notes

`PolynomialQ[expr, var]` tests whether `expr` is a polynomial in `var` — i.e.
expands to a sum of products of non-negative integer powers of the variables
with variable-free coefficients. `Sin[x] + x` fails because of the
transcendental term, and `x^2 + y/x` fails because `y/x = y x^(-1)` carries a
negative power of `x`. The multivariate form `PolynomialQ[expr, {x, y}]`
requires polynomiality in all listed variables simultaneously. Note that the
test is syntactic up to expansion and does not cancel rational expressions: a
fraction such as `(x^2 - 1)/(x - 1)` is reported `False` even though it equals
the polynomial `x + 1`.
