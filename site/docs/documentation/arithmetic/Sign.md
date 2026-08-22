# Sign

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Sign[x] gives -1, 0, or 1 for real numeric x according to its sign, and z/Abs[z] for a nonzero numeric complex z.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (6)

```mathematica
In[1]:= Sign[-7]
Out[1]= -1

In[2]:= Sign[0]
Out[2]= 0

In[3]:= Sign[{-2, 0, 5}]
Out[3]= {-1, 0, 1}

In[4]:= Sign[3 + 4 I]
Out[4]= 3/5 + 4/5*I

In[5]:= Sign[(1 + I)^2]
Out[5]= I

In[6]:= Sign[2 - 2 I]
Out[6]= (1/2 - 1/2*I) Sqrt[2]
```

## Implementation notes

`builtin_sign` returns the sign (-1/0/1) of a real number — direct comparisons for `EXPR_INTEGER`/`EXPR_REAL` (and Rational by sign of numerator×denominator), `mpz_sgn` for BigInt, `mpfr_sgn` for MPFR. For a numeric `Complex[re, im]` with both parts numeric it returns the unit-modulus direction `z/Abs[z]` (short-circuiting `0+0I -> 0`); MPFR components take a fast path computing the direction directly via `mpfr_hypot` and division at the combined working precision rather than building the symbolic `z·Power[Abs[z], -1]` tree. Non-numeric arguments return `NULL` (unevaluated).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [ReIm](../../arithmetic/ReIm/), [Abs](../../arithmetic/Abs/), [Conjugate](../../arithmetic/Conjugate/), [Arg](../../arithmetic/Arg/), [Rational](../../arithmetic/Rational/), [Complex](../../arithmetic/Complex/)

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_bignum_rational_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bignum_rational_numeric.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)

## Notes & additional examples

### Notes

For real `x`, `Sign[x]` is -1, 0, or 1. For a nonzero complex `z` it returns the unit-modulus direction `z/Abs[z]`: `(1+I)^2 = 2I` points straight up the imaginary axis (`I`), while `2 - 2I` lies on the diagonal and returns the exact unit vector `(1 - I)/Sqrt[2]`. Sign is Listable.
