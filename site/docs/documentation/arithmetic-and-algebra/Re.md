# Re

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Re[z] gives the real part of numeric z; Re[Re[z]], Re[Im[z]], Re[Abs[z]], Re[Arg[z]] fold since those heads are real-valued.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_re` returns the real part. It returns the argument itself for real numeric kinds (`EXPR_INTEGER`/`EXPR_REAL`/`EXPR_MPFR`/Rational) and for real-valued head calls (`Re`/`Im`/`Abs`/`Arg`, via `is_real_valued_head_call`); for a `Complex[re, im]` literal it returns `re`; and for an expression `complex_decompose` splits into numeric real/imaginary parts, it returns the real part. Otherwise (genuinely symbolic) it returns `NULL` and the call stays unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
