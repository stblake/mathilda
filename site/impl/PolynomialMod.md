---
source: src/poly/poly.c
---
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
