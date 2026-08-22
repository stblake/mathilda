# Arg

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Arg[z] gives the argument (phase angle in (-Pi, Pi]) of numeric z; 0 for nonnegative reals, Pi for negative reals.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (7)

```mathematica
In[1]:= Arg[1]
Out[1]= 0

In[2]:= Arg[-1]
Out[2]= Pi

In[3]:= Arg[1 + I]
Out[3]= 1/4 Pi

In[4]:= Arg[(1 + I)^10]
Out[4]= 1/2 Pi

In[5]:= Arg[-2 + 2 I]
Out[5]= 3/4 Pi

In[6]:= Table[Arg[(1 + I)^k], {k, 0, 8}]
Out[6]= {0, 1/4 Pi, 1/2 Pi, 3/4 Pi, Pi, -3/4 Pi, -1/2 Pi, -1/4 Pi, 0}

In[7]:= N[Arg[2 + 3 I], 40]
Out[7]= 0.98279372324732906798571061101466601449686
```

## Implementation notes

`builtin_arg` returns the phase angle in `(-Pi, Pi]`. A pure MPFR real folds to exact `0` or `Pi` by sign. For a `Complex[re, im]` whose parts are exact (Integer/Rational), it recognises the special directions and returns exact multiples of `Pi`: `0` for positive reals, `Pi` for negatives, `±Pi/2` on the imaginary axis, and `±Pi/4`, `±3Pi/4` on the diagonals; otherwise it returns the symbolic `ArcTan[re, im]`. When either component carries MPFR it evaluates `mpfr_atan2` at the combined precision; an inexact machine `Real` falls through to the libm `atan2(im, re)`. Symbolic inputs return `NULL`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [ReIm](../../arithmetic/ReIm/), [Abs](../../arithmetic/Abs/), [Sign](../../arithmetic/Sign/), [Conjugate](../../arithmetic/Conjugate/), [Rational](../../arithmetic/Rational/), [Complex](../../arithmetic/Complex/)

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_bignum_rational_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bignum_rational_numeric.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)

## Notes & additional examples

### Notes

`Arg[z]` gives the phase angle in the range `(-Pi, Pi]`: 0 for positive reals, `Pi` for negative reals.
