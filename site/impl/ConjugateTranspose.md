---
source: src/list.c
---
`builtin_conjugate_transpose` (in `src/list.c`) is a thin composition over existing primitives. It first checks the argument is a rectangular nested `List` via `get_array_dimensions`; a symbolic (non-list) matrix is left unevaluated so `ConjugateTranspose[A]` survives. For a 1-D vector it just maps `Conjugate` elementwise. Otherwise it builds and evaluates `Transpose[m]` (or `Transpose[m, spec]`), then conjugates the transposed result. All heavy lifting is delegated to `Transpose` and `Conjugate` through `eval_and_free`.
