---
source: src/complex.c
---
`builtin_reim` returns `{Re[z], Im[z]}` as a two-element `List`. It mirrors `builtin_re`/`builtin_im`: `{f[z], 0}` for real-valued head calls (`Re`/`Im`/`Abs`/`Arg`), `{re, im}` for a `Complex[re, im]` literal, `{x, 0}` for real numeric kinds (Integer/Real/Rational), and `{re, im}` when `complex_decompose` yields numeric parts. Returns `NULL` (unevaluated) for genuinely symbolic input.
