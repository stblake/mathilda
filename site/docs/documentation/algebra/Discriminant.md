# Discriminant

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Discriminant[poly, var]
    gives the discriminant of poly with respect to var: up to sign and
    leading-coefficient scaling, Resultant[poly, D[poly, var], var] /
    lc[poly, var].  Vanishes iff poly has a repeated root in var.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Discriminant[a x^2 + b x + c, x]
Out[1]= b^2 - 4 a c

In[2]:= Discriminant[5 x^4 - 3 x + 9, x]
Out[2]= 23273325

In[3]:= Discriminant[(x-1)(x-2)(x-3), x]
Out[3]= 4

In[4]:= Discriminant[(x-1)(x-2)(x-1), x]
Out[4]= 0
```

## Implementation notes

**Algorithm.** `builtin_discriminant` (in `src/poly/poly.c`) computes `Disc_x(p)` from the standard resultant identity. After confirming `p` is a polynomial in `x` (`PolynomialQ`) and expanding it, it reads the degree `n` (returning 0 for degree 0/1). It then forms the derivative `p'` with `poly_derivative`, computes `R = Res_x(p, p')` via `resultant_internal`, and applies the closed form

  Disc = (-1)^{n(n-1)/2} · Res(p, p') / a_n,

where `a_n` is the leading coefficient (obtained with `get_coeff`). The quotient is built as `Times[sign·R, a_n^{-1}]` and the final result is `expr_expand`ed.

**Data structures.** Plain `Expr*` polynomials throughout; the heavy lifting is the subresultant/Euclidean `resultant_internal` (see `Resultant`). Degree and coefficient queries use `get_degree_poly`/`get_coeff`.

**Complexity / limits.** Dominated by the resultant computation, which is polynomial in `n` and the coefficient sizes but can be expensive for large multivariate inputs.

- `Protected`, `Listable`.
- Computes the discriminant of polynomial `poly` with respect to `var`.
- The discriminant is zero if and only if the polynomial has multiple roots.
- Derived symbolically utilizing the formula $D = \frac{(-1)^{n(n-1)/2}}{a_n} Resultant(P, P', var)$.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- B. L. van der Waerden, *Algebra*, vol. 1 (Springer).
- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Discriminant[x^2 - 5 x + 6, x]
Out[1]= 1
```

```mathematica
In[1]:= Discriminant[a x^2 + b x + c, x]
Out[1]= b^2 - 4 a c
```

```mathematica
In[1]:= Discriminant[x^3 + p x + q, x]
Out[1]= -4 p^3 - 27 q^2
```

```mathematica
In[1]:= Discriminant[x^4 + 1, x]
Out[1]= 256
```

### Notes

`Discriminant[poly, var]` is, up to sign and leading-coefficient scaling,
`Resultant[poly, D[poly, var], var] / lc[poly, var]`, and it vanishes exactly when
`poly` has a repeated root in `var`. The first example has distinct roots `2` and
`3`, so its discriminant is a nonzero constant. The second recovers the familiar
`b^2 - 4 a c` from the quadratic formula. The third gives the classical
depressed-cubic discriminant `-4 p^3 - 27 q^2`, whose sign distinguishes three real
roots from one real and two complex roots. The fourth shows `x^4 + 1`, which has
four distinct complex roots and discriminant `256`.
