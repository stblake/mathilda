---
source: src/poly/algebraicnumberpolynomial.c
references:
  - "H. Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138 (Springer, 1993), §4.2 (power-basis representation of number-field elements)."
---
**Algorithm.** `AlgebraicNumberPolynomial[a, x]` is a purely **structural** read
and rebuild — no field arithmetic and no FLINT. The coefficient vector of an
`AlgebraicNumber` object already lives in the object, so the defining polynomial
is assembled directly:

1. An `Integer`, `BigInt`, or `Rational` `a` is the constant polynomial and is
   returned unchanged (`expr_copy`).
2. For `a = AlgebraicNumber[theta, {c0, c1, …, cn}]` whose coefficients are all
   exact rationals (`Integer`/`BigInt`/`Rational` — the forms canonicalisation
   stores), it builds `c0 + c1 x + c2 x^2 + … + cn x^n` as a `Plus` of `Times`/
   `Power` terms. The generator `theta` is irrelevant here: `a` is recovered by
   substituting `x -> theta`. An empty coefficient list yields `0`.
3. Any other argument routes an `AlgebraicNumberPolynomial::naobj` message
   through `mth_message` and stays unevaluated.

The returned `Plus`/`Times`/`Power` tree is left **unevaluated**; the evaluator
canonicalises it on its next fixed-point step (zero terms drop, `Times[1, x]`
folds to `x`, monomials sort by degree).

**Data structures.** `Expr` trees only. It reads the coefficient `List` out of
the `AlgebraicNumber[theta, {coeffs}]` representation and emits a polynomial
`Expr`; the `anp_is_rational_coeff` guard enforces the exact-rational invariant
of a canonical object.

**Complexity / limits.** `O(n)` in the degree — one monomial per coefficient,
with a single allocation sweep. `Listable`, `Protected`. It declines on anything
that is not a rational or a well-formed `AlgebraicNumber` object (including an
`AlgebraicNumber` with a symbolic or inexact coefficient).
