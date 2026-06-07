# Im

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Im[z] gives the imaginary part of numeric z, and 0 for real or real-valued (Re/Im/Abs/Arg) arguments.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_im` returns the imaginary part. It returns `0` for the real-valued-by-construction heads (`Re`/`Im`/`Abs`/`Arg`) and for any real numeric kind (Integer/Real/Rational/MPFR), copies the second component of a `Complex[re, im]` literal, and for a general expression runs `complex_decompose` (a recursive Plus/Times walk that propagates `Complex` literals through complex multiplication) — returning the imaginary part only when both decomposed parts are concretely numeric (`is_numeric_real`). Otherwise `NULL`, leaving the symbolic head in place. `Re`/`ReIm` in the same file share this machinery.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Im[3 + 4 I]
Out[1]= 4

In[2]:= Im[7]
Out[2]= 0
```

### Notes

`Im[z]` extracts the imaginary part of numeric `z`, giving 0 for real (or real-valued) arguments. It is Listable.
