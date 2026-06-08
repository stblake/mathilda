# Sign

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Sign[x] gives -1, 0, or 1 for real numeric x according to its sign, and z/Abs[z] for a nonzero numeric complex z.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_sign` returns the sign (-1/0/1) of a real number — direct comparisons for `EXPR_INTEGER`/`EXPR_REAL` (and Rational by sign of numerator×denominator), `mpz_sgn` for BigInt, `mpfr_sgn` for MPFR. For a numeric `Complex[re, im]` with both parts numeric it returns the unit-modulus direction `z/Abs[z]` (short-circuiting `0+0I -> 0`); MPFR components take a fast path computing the direction directly via `mpfr_hypot` and division at the combined working precision rather than building the symbolic `z·Power[Abs[z], -1]` tree. Non-numeric arguments return `NULL` (unevaluated).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Sign[-7]
Out[1]= -1

In[2]:= Sign[0]
Out[2]= 0

In[3]:= Sign[3 + 4 I]
Out[3]= 3/5 + 4/5*I

In[4]:= Sign[{-2, 0, 5}]
Out[4]= {-1, 0, 1}
```

### Notes

For real `x`, `Sign[x]` is -1, 0, or 1. For a nonzero complex `z` it returns the unit-modulus direction `z/Abs[z]`. Sign is Listable.
