---
source: src/arithmetic.c
---
`builtin_complex` (`src/arithmetic.c`) is the `Complex[re, im]` constructor's auto-simplifier. It only collapses zero-imaginary cases: `Complex[r, 0]` (integer 0) returns `re` unchanged; `Complex[r, 0.0]` (real 0) returns `re`, promoting an integer real part to a `Real` so the result stays inexact. Any genuinely complex value returns `NULL`, leaving the literal `Complex[re, im]` in place as the canonical representation. Re/Im decomposition, arithmetic on complex values, and printing as `a + b I` live elsewhere (`src/complex.c`, `src/print.c`); this handler is purely the constructor normalisation step.
