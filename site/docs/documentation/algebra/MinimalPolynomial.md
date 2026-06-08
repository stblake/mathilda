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

_No verified examples yet for this function._

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
