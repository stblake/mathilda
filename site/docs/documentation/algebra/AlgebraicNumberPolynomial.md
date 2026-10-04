# AlgebraicNumberPolynomial

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AlgebraicNumberPolynomial[a, x]`**

gives the polynomial in x corresponding to the AlgebraicNumber object a.

<details>
<summary>Notes</summary>

For a = AlgebraicNumber\[theta, {c0, c1, ..., cn}\], the result is the polynomial c0 + c1 x + ... + cn x^n, from which a is recovered by replacing x with theta.  An integer or rational a is the constant polynomial and is returned unchanged; any other argument stays unevaluated.  Threads over lists.  Attributes: Listable, Protected.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2] + Sqrt[3], {1, 2, 3, 4}], x]
Out[1]= 1 + 2 x + 3 x^2 + 4 x^3
```

Threads

```mathematica
In[2]:= AlgebraicNumberPolynomial[{2, AlgebraicNumber[Sqrt[2], {1, 2}]}, x]
Out[2]= {2, 1 + 2 x}
```

```mathematica
In[3]:= a = AlgebraicNumber[Sqrt[2 + Sqrt[3]], {1, 2, 3, 4}]; RootReduce[(AlgebraicNumberPolynomial[a, x] /. x -> Sqrt[2 + Sqrt[3]]) == a]
Out[3]= True
```

### Worked examples (2)

```mathematica
In[4]:= AlgebraicNumberPolynomial[2, x]
Out[4]= 2

In[5]:= AlgebraicNumberPolynomial[1/2, x]
Out[5]= 1/2
```

### Applications (2)

C0 + c1 x from {1, 2}

```mathematica
In[6]:= AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2], {1, 2}], x]
Out[6]= 1 + 2 x
```

A rational is the constant polynomial

```mathematica
In[7]:= AlgebraicNumberPolynomial[7, x]
Out[7]= 7
```

## Algorithm

algebraicnumberpolynomial.c — AlgebraicNumberPolynomial[a, x].

See algebraicnumberpolynomial.h for the contract. The coefficient vector of an AlgebraicNumber object is stored in the object itself, so producing the defining polynomial c0 + c1 x + ... + cn x^n is a purely structural read + build: no field arithmetic and no FLINT. The returned Plus/Times/Power tree is canonicalised by the evaluator on its next fixed-point step (zero terms drop, Times[1, x] folds to x, monomials sort by degree).

## Implementation notes

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

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [AlgebraicNumber](../../algebra/AlgebraicNumber/), [RootReduce](../../algebra/RootReduce/)

- H. Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138 (Springer, 1993), §4.2 (power-basis representation of number-field elements).
- Source: [`src/poly/algebraicnumberpolynomial.c`](https://github.com/stblake/mathilda/blob/main/src/poly/algebraicnumberpolynomial.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_algebraicnumberpolynomial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumberpolynomial.c)

## Notes & additional examples

### Notes

`AlgebraicNumberPolynomial[a, x]` gives the polynomial in `x` whose coefficients
are the stored coefficient list of the `AlgebraicNumber` object `a`: for
`a = AlgebraicNumber[theta, {c0, c1, …, cn}]` the result is
`c0 + c1 x + … + cn x^n`, and `a` is recovered by replacing `x` with `theta`. The
generator `theta` plays no part — this is a purely structural read of the
coefficient vector, with no field arithmetic.

An integer or rational `a` is the constant polynomial and is returned unchanged;
any other argument stays unevaluated. `Listable` and `Protected`.
