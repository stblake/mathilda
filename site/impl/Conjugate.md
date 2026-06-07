---
source: src/complex.c
---
`builtin_conjugate` folds the involution `Conjugate[Conjugate[z]] -> z` and treats the real-valued-by-construction heads `Re`, `Im`, `Abs`, `Arg` as fixed points. For a `Complex[re, im]` literal it returns `make_complex(re, -im)`; for real numerics (Integer/Real/Rational) and any expression that `is_numeric_real` (e.g. `Sqrt[2]`, `Pi`) it returns the argument unchanged. An expression that `complex_decompose` splits into concretely-numeric real/imag parts is conjugated as `re - im*I`. Symbolic inputs return `NULL` (a one-argument arity check emits `Conjugate::argx`).
