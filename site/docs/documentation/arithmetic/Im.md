# Im

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Im[z] gives the imaginary part of numeric z, and 0 for real or real-valued (Re/Im/Abs/Arg) arguments.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (8)

```mathematica
In[1]:= Im[3 + 4 I]
Out[1]= 4

In[2]:= Im[7]
Out[2]= 0

In[3]:= Im[Sqrt[-4]]
Out[3]= 2

In[4]:= Im[(1 + I)^10]
Out[4]= 32

In[5]:= Im[Log[-1]]
Out[5]= Pi

In[6]:= Im[{1 + 2 I, 3 - 4 I, 5}]
Out[6]= {2, -4, 0}

In[7]:= Im[Gamma[1 + I]]
Out[7]= Im[Gamma[1 + I]]

In[8]:= N[Im[Gamma[1 + I]], 40]
Out[8]= -0.15494982830181068512495513048388660519589
```

## Implementation notes

`builtin_im` returns the imaginary part. It returns `0` for the real-valued-by-construction heads (`Re`/`Im`/`Abs`/`Arg`) and for any real numeric kind (Integer/Real/Rational/MPFR), copies the second component of a `Complex[re, im]` literal, and for a general expression runs `complex_decompose` (a recursive Plus/Times walk that propagates `Complex` literals through complex multiplication) — returning the imaginary part only when both decomposed parts are concretely numeric (`is_numeric_real`). Otherwise `NULL`, leaving the symbolic head in place. `Re`/`ReIm` in the same file share this machinery.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Re](../../arithmetic/Re/), [ReIm](../../arithmetic/ReIm/), [Abs](../../arithmetic/Abs/), [Sign](../../arithmetic/Sign/), [Conjugate](../../arithmetic/Conjugate/), [Arg](../../arithmetic/Arg/), [Rational](../../arithmetic/Rational/), [Complex](../../arithmetic/Complex/)

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_airyai.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airyai.c)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)

## Notes & additional examples

### Notes

`Im[z]` extracts the imaginary part of numeric `z`, giving 0 for real (or real-valued) arguments. It is Listable.
