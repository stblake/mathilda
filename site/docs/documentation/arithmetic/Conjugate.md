# Conjugate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Conjugate[z] gives the complex conjugate Re[z] - I Im[z] of numeric z; real and real-valued (Re/Im/Abs/Arg) arguments are returned unchanged.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= Conjugate[3 + 4 I]
Out[1]= 3 - 4*I

In[2]:= Conjugate[5]
Out[2]= 5

In[3]:= Conjugate[{1 + I, 2 - 3 I}]
Out[3]= {1 - I, 2 + 3*I}

In[4]:= Conjugate[(2 + I)/(1 - 3 I)]
Out[4]= -1/10 - 7/10*I

In[5]:= z Conjugate[z] /. z -> 3 + 4 I
Out[5]= 25
```

## Implementation notes

`builtin_conjugate` folds the involution `Conjugate[Conjugate[z]] -> z` and treats the real-valued-by-construction heads `Re`, `Im`, `Abs`, `Arg` as fixed points. For a `Complex[re, im]` literal it returns `make_complex(re, -im)`; for real numerics (Integer/Real/Rational) and any expression that `is_numeric_real` (e.g. `Sqrt[2]`, `Pi`) it returns the argument unchanged. An expression that `complex_decompose` splits into concretely-numeric real/imag parts is conjugated as `re - im*I`. Symbolic inputs return `NULL` (a one-argument arity check emits `Conjugate::argx`).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [ReIm](../../arithmetic/ReIm/), [Abs](../../arithmetic/Abs/), [Sign](../../arithmetic/Sign/), [Arg](../../arithmetic/Arg/), [Rational](../../arithmetic/Rational/), [Complex](../../arithmetic/Complex/)

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_blas.c`](https://github.com/stblake/mathilda/blob/main/tests/test_blas.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_complexexpand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_complexexpand.c)
- Tests: [`tests/test_conjugate_transpose.c`](https://github.com/stblake/mathilda/blob/main/tests/test_conjugate_transpose.c)

## Notes & additional examples

### Notes

`Conjugate[z]` returns `Re[z] - I Im[z]`; real arguments are returned unchanged. It is Listable.
