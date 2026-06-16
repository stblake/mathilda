# MinimalPolynomial

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
MinimalPolynomial[s, x]
    gives the lowest-degree polynomial in x with integer coefficients,
    positive leading coefficient and content 1, having the algebraic
    number s as a root.  s may be built from rationals, radicals, the
    imaginary unit, roots of unity, and Root[] objects.
MinimalPolynomial[s]
    gives the minimal polynomial as a pure function.
MinimalPolynomial[s, x, Extension -> a]
    gives the characteristic polynomial of s in Q(a) over Q(a).
    Computed by resultant elimination of the radicals; threads over
    lists.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= MinimalPolynomial[Sqrt[2] + Sqrt[3], x]
Out[1]= 1 - 10 x^2 + x^4

In[2]:= MinimalPolynomial[(1 + I)/Sqrt[2], x]
Out[2]= 1 + x^4

In[3]:= MinimalPolynomial[Root[2 #1^3 - 2 #1 + 7 &, 1] + 17, x]
Out[3]= -9785 + 1732 x - 102 x^2 + 2 x^3

In[4]:= MinimalPolynomial[Sqrt[2], x, Extension -> E^(I Pi/4)]
Out[4]= 4 - 4 x^2 + x^4
```

## Implementation notes

**Algorithm.** `builtin_minimalpolynomial` computes the minimal polynomial of an algebraic number `s` over Q in `x` via resultant elimination:

1. *Atom walk* (`mp_walk`): recursively rewrite `s` into a "value expression" V in fresh auxiliary symbols `t_i`, recording for each `t_i` a polynomial defining relation `p_i` (in `t_i` and earlier auxiliaries). Every relation is kept polynomial — negative powers are turned into reciprocal variables (`D·w - 1`) so no fractions appear.
2. *Elimination*: build `g = (x - V)` and eliminate each `t_i` (highest index first) by `g <- Resultant[g, p_i, t_i]`. The introduction order guarantees `t_i` is present in `g` when eliminated and that each relation references only earlier auxiliaries, so the chain terminates in a univariate `G(x)`.
3. *Clear and factor*: take `Numerator[Together[G]]`, make primitive over Z, and factor. Evaluate `s` numerically to high precision and select the unique irreducible factor that vanishes at `s`.
4. Return that factor, primitive with positive leading coefficient.

`Extension -> a` uses the tower law: if `s ∈ Q(a)`, the characteristic polynomial of `s` over Q(a) is `m_s(x)^([Q(a):Q]/[Q(s):Q])`; membership is checked via the primitive-element degree `[Q(a,s):Q] == [Q(a):Q]`.

**Data structures.** `Expr*` trees with fresh internal auxiliary symbols; resultants run on the multivariate polynomial machinery (`Resultant`), denominator clearing through `internal_together`/`Numerator`, primitivisation/factoring via the Z-polynomial routines (`zupoly`, `facpoly`), and root selection through high-precision `numericalize` (MPFR when built). The builtin takes ownership of `res` and returns a fresh `Expr*` or NULL.

**Complexity / limits.** Dominated by the iterated resultant elimination (one resultant per auxiliary symbol, with the usual degree blow-up) and the final univariate factorisation. Numeric root-matching disambiguates the irreducible factor.

- `Listable`, `Protected`. A `List` first argument threads element-wise.
- `s` may be built from integers and rationals, radicals (`Sqrt`,

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/poly/minpoly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/minpoly.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= MinimalPolynomial[Sqrt[2], x]
Out[1]= -2 + x^2
```

```mathematica
In[1]:= MinimalPolynomial[(1 + Sqrt[5])/2, x]
Out[1]= -1 - x + x^2
```

```mathematica
In[1]:= MinimalPolynomial[Sqrt[2] + Sqrt[3], x]
Out[1]= 1 - 10 x^2 + x^4
```

```mathematica
In[1]:= MinimalPolynomial[Cos[2 Pi/5], x]
Out[1]= -1 + 2 x + 4 x^2
```

```mathematica
In[1]:= MinimalPolynomial[Sqrt[2 + Sqrt[2]], x]
Out[1]= 2 - 4 x^2 + x^4
```

### Notes

`MinimalPolynomial[s, x]` returns the lowest-degree integer polynomial in `x`,
with positive leading coefficient and content 1, that has the algebraic number
`s` as a root. The first two examples recover the defining polynomials of `Sqrt[2]`
and the golden ratio (`x^2 - x - 1`). The real power shows up with compound
algebraic numbers: `Sqrt[2] + Sqrt[3]` is degree 4 over the rationals, and
`MinimalPolynomial` finds its quartic `x^4 - 10 x^2 + 1` by eliminating the
radicals with resultants — not by numerical root-finding. It also handles
algebraic constants beyond plain radicals: `Cos[2 Pi/5]` (a root of a cyclotomic
relation) yields `4 x^2 + 2 x - 1`, and the nested radical `Sqrt[2 + Sqrt[2]]`
yields the quartic `x^4 - 4 x^2 + 2`. The input may be built from rationals,
radicals, the imaginary unit, roots of unity, and `Root[]` objects.
