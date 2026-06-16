# Decompose

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Decompose[poly, x]
    decomposes the univariate polynomial poly into the deepest possible
    composition {p1, p2, ..., pk} such that poly == p1[p2[...[pk[x]]...]],
    with each pi a polynomial of degree >= 2 in x.
    Returns {poly} if no nontrivial decomposition exists.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_decompose` (`src/poly/poly.c`) finds polynomials with f(x) = g(h(x)), deg(h) ≥ 2. It first verifies the input is a polynomial in x via `internal_polynomialq` (returns `NULL`/unevaluated otherwise) and then calls the recursive driver `decompose_recursive`. That driver expands f, takes its degree n, and applies two reductions in order: (1) a common-degree shortcut — if the gcd `d` of the degrees of all nonzero monomials exceeds 1, every term is divisible by x^d, so it substitutes y = x^d and recurses; (2) trial composition (the Kozen–Landau scheme) — for each proper divisor s of n it builds the candidate inner polynomial of degree s by matching the high-order coefficients of `a_n*(x^s)^r` against f (solving for each coefficient with `exact_poly_div` over the non-x variables), then divides f by successive powers of the candidate inner h with `poly_div_rem`, accepting the decomposition only if every remainder is x-free and the final quotient is zero. The outer factor list `Lg` from the recursion is concatenated with the inner h to form the returned `List`. If no decomposition is found it returns `{f}`.

**Data structures.** `Expr*` polynomial trees throughout; coefficients are extracted with `get_coeff`/`get_all_coeffs_expanded` and variables collected via `collect_variables`. Exact divisions use the multivariate `exact_poly_div` and `poly_div_rem` helpers in the same file.

**Complexity / limits.** Univariate decomposition over the polynomial's coefficient domain; multivariate coefficients are tolerated through the variable-list machinery. Cost is dominated by repeated `expr_expand` and `poly_div_rem` over the divisors of n.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- D. Kozen and S. Landau, "Polynomial Decomposition Algorithms", J. Symbolic Computation 7 (1989).
- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Decompose[x^4 + 6 x^2 + 10, x]
Out[1]= {10 + 6 x + x^2, x^2}
```

```mathematica
In[1]:= Decompose[x^6 + 6 x^4 + 12 x^2 + 8, x]
Out[1]= {8 + 12 x + 6 x^2 + x^3, x^2}
```

```mathematica
In[1]:= Decompose[x^4 + 2 x^3 + 3 x^2 + 2 x + 5, x]
Out[1]= {5 + 2 x + x^2, x + x^2}
```

```mathematica
In[1]:= Decompose[x^6, x]
Out[1]= {x^6}
```

### Notes

`Decompose[poly, x]` returns `{p1, p2, ..., pk}` such that `poly == p1[p2[...[pk[x]]...]]`, with each `pi` of degree >= 2. Out[1] reads as `p1 = 10 + 6 x + x^2` composed with `p2 = x^2`. The inner factor need not be a pure power: in the third example `Decompose` recognises the non-monomial inner map `x + x^2`, recovering `(x^2 + x)^2 + 2 (x^2 + x) + 5`. A polynomial with no nontrivial decomposition (for example a single monomial like `x^6`) returns `{poly}` unchanged.
