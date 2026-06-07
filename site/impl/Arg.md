---
source: src/complex.c
---
`builtin_arg` returns the phase angle in `(-Pi, Pi]`. A pure MPFR real folds to exact `0` or `Pi` by sign. For a `Complex[re, im]` whose parts are exact (Integer/Rational), it recognises the special directions and returns exact multiples of `Pi`: `0` for positive reals, `Pi` for negatives, `±Pi/2` on the imaginary axis, and `±Pi/4`, `±3Pi/4` on the diagonals; otherwise it returns the symbolic `ArcTan[re, im]`. When either component carries MPFR it evaluates `mpfr_atan2` at the combined precision; an inexact machine `Real` falls through to the libm `atan2(im, re)`. Symbolic inputs return `NULL`.
