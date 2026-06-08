# PolynomialMod

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PolynomialMod[poly, m]
    reduces poly modulo m.  If m is an integer, each coefficient of
    poly is reduced to a canonical residue in {0, ..., m-1}.  If m is a
    polynomial, poly is reduced modulo m as polynomials over the
    rationals (in contrast to PolynomialRemainder, the leading
    coefficient of m is not normalised).
PolynomialMod[poly, {m1, m2, ...}]
    reduces modulo each mi in turn.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PolynomialMod[3x^2+2x+1,2]
Out[1]= 1 + x^2

In[2]:= PolynomialMod[3x^2+2x+1,x^2+1]
Out[2]= -2 + 2 x

In[3]:= PolynomialMod[35x^3+21x^2 y^2-17x y^3+55z-123,19]
Out[3]= 10 + 16 x^3 + 2 x^2 y^2 + 2 x y^3 + 17 z

In[4]:= PolynomialMod[3x^3+21x^2 y^2-7x y^3+55,{2x^2-7,x y-3, 9}]
Out[4]= 1 + 7 x + x^3 + 4 y^2
```

## Implementation notes

**Algorithm.** `builtin_polynomialmod` reduces a polynomial modulo `m`, where `m` may be an
integer (reduce each coefficient mod m, into the symmetric/least-residue range), a polynomial
(polynomial remainder), or a `List` of moduli applied successively. It threads over structural
heads (`List`, `Equal`, `Less`, `And`, `Not`, …) by recursing into each argument. The actual
reduction is done by `polynomial_mod_single`; when `m` is a list with an integer member the
integer reduction is applied alongside the polynomial reductions. Unlike `PolynomialRemainder`,
`PolynomialMod` performs no leading-coefficient division — coefficients are reduced rather than
divided — so it stays within the coefficient ring (the modular-arithmetic analogue of `Mod` on
integers, lifted coefficientwise).

**Data structures.** Ordinary `Expr` polynomial trees; integer coefficient reduction uses the
standard integer/bigint arithmetic helpers.

- `Protected`, `Listable`.
- Reduces a polynomial modulo an integer, another polynomial, or a list of integers/polynomials.
- Always gives a result with minimal degree and leading coefficients.
- Handles rational division mapping perfectly scaling over exact modulo structures dynamically.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
