# Complex

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Complex[re, im]`**

represents the complex number re + im I.

<details>
<summary>Notes</summary>

Complex is the canonical head produced by arithmetic when an Integer, Real, or Rational acquires an imaginary part. Pure-real inputs collapse to the underlying number; im == 0 unwraps to re.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= Complex[3, 4]
Out[1]= 3 + 4*I

In[2]:= Complex[5, 0]
Out[2]= 5

In[3]:= (1 + I)^8
Out[3]= 16

In[4]:= (2 + 3 I)/(1 - I)
Out[4]= -1/2 + 5/2*I

In[5]:= N[(3 + 4 I)^(1/3), 40]
Out[5]= 1.6289371459221758752146093717175049715341 + 0.52017450230454583954569417015944746788805*I
```

## Implementation notes

`builtin_complex` (`src/arithmetic.c`) is the `Complex[re, im]` constructor's auto-simplifier. It only collapses zero-imaginary cases: `Complex[r, 0]` (integer 0) returns `re` unchanged; `Complex[r, 0.0]` (real 0) returns `re`, promoting an integer real part to a `Real` so the result stays inexact. Any genuinely complex value returns `NULL`, leaving the literal `Complex[re, im]` in place as the canonical representation. Re/Im decomposition, arithmetic on complex values, and printing as `a + b I` live elsewhere (`src/complex.c`, `src/print.c`); this handler is purely the constructor normalisation step.

**Attributes:** none registered.

## References

- Source: [`src/arithmetic.c`](https://github.com/stblake/mathilda/blob/main/src/arithmetic.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_beta.c`](https://github.com/stblake/mathilda/blob/main/tests/test_beta.c)
- Tests: [`tests/test_blas.c`](https://github.com/stblake/mathilda/blob/main/tests/test_blas.c)
- Tests: [`tests/test_chop.c`](https://github.com/stblake/mathilda/blob/main/tests/test_chop.c)

## Notes & additional examples

### Notes

`Complex[re, im]` is the canonical head for `re + im I`; a zero imaginary part collapses back to the underlying real number.
