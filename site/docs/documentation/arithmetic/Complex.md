# Complex

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Complex[re, im]
    represents the complex number re + im I.
Complex is the canonical head produced by arithmetic when an Integer,
Real, or Rational acquires an imaginary part. Pure-real inputs collapse
to the underlying number; im == 0 unwraps to re.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_complex` (`src/arithmetic.c`) is the `Complex[re, im]` constructor's auto-simplifier. It only collapses zero-imaginary cases: `Complex[r, 0]` (integer 0) returns `re` unchanged; `Complex[r, 0.0]` (real 0) returns `re`, promoting an integer real part to a `Real` so the result stays inexact. Any genuinely complex value returns `NULL`, leaving the literal `Complex[re, im]` in place as the canonical representation. Re/Im decomposition, arithmetic on complex values, and printing as `a + b I` live elsewhere (`src/complex.c`, `src/print.c`); this handler is purely the constructor normalisation step.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/arithmetic.c`](https://github.com/stblake/mathilda/blob/main/src/arithmetic.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Complex[3, 4]
Out[1]= 3 + 4*I

In[2]:= (3 + 4 I) + (1 - 2 I)
Out[2]= 4 + 2*I

In[3]:= Complex[5, 0]
Out[3]= 5
```

### Notes

`Complex[re, im]` is the canonical head for `re + im I`; a zero imaginary part collapses back to the underlying real number.
