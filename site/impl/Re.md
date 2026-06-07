---
source: src/complex.c
---
`builtin_re` returns the real part. It returns the argument itself for real numeric kinds (`EXPR_INTEGER`/`EXPR_REAL`/`EXPR_MPFR`/Rational) and for real-valued head calls (`Re`/`Im`/`Abs`/`Arg`, via `is_real_valued_head_call`); for a `Complex[re, im]` literal it returns `re`; and for an expression `complex_decompose` splits into numeric real/imaginary parts, it returns the real part. Otherwise (genuinely symbolic) it returns `NULL` and the call stays unevaluated.
