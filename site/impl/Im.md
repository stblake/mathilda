---
source: src/complex.c
---
`builtin_im` returns the imaginary part. It returns `0` for the real-valued-by-construction heads (`Re`/`Im`/`Abs`/`Arg`) and for any real numeric kind (Integer/Real/Rational/MPFR), copies the second component of a `Complex[re, im]` literal, and for a general expression runs `complex_decompose` (a recursive Plus/Times walk that propagates `Complex` literals through complex multiplication) — returning the imaginary part only when both decomposed parts are concretely numeric (`is_numeric_real`). Otherwise `NULL`, leaving the symbolic head in place. `Re`/`ReIm` in the same file share this machinery.
