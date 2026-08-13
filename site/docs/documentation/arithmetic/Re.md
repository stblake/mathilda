# Re

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Re[z] gives the real part of numeric z; Re[Re[z]], Re[Im[z]], Re[Abs[z]], Re[Arg[z]] fold since those heads are real-valued.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= Re[3 + 4 I]
Out[1]= 3

In[2]:= Re[7]
Out[2]= 7

In[3]:= Re[{1 + I, 2 - 3 I, 5}]
Out[3]= {1, 2, 5}

In[4]:= Re[(1 + I)^10]
Out[4]= 0

In[5]:= Re[1/(2 + 3 I)]
Out[5]= 2/13
```

## Implementation notes

`builtin_re` returns the real part. It returns the argument itself for real numeric kinds (`EXPR_INTEGER`/`EXPR_REAL`/`EXPR_MPFR`/Rational) and for real-valued head calls (`Re`/`Im`/`Abs`/`Arg`, via `is_real_valued_head_call`); for a `Complex[re, im]` literal it returns `re`; and for an expression `complex_decompose` splits into numeric real/imaginary parts, it returns the real part. Otherwise (genuinely symbolic) it returns `NULL` and the call stays unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Im](../../arithmetic/Im/), [ReIm](../../arithmetic/ReIm/), [Abs](../../arithmetic/Abs/), [Sign](../../arithmetic/Sign/), [Conjugate](../../arithmetic/Conjugate/), [Arg](../../arithmetic/Arg/), [Rational](../../arithmetic/Rational/), [Complex](../../arithmetic/Complex/)

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_airyai.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airyai.c)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)

## Notes & additional examples

### Notes

`Re[z]` extracts the real part of numeric `z`; a purely real argument is returned unchanged. It is Listable.
