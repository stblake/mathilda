---
source: src/power.c
---
`builtin_sqrt` is a thin wrapper: it rewrites `Sqrt[x]` to `Power[x, Rational[1, 2]]` (via `make_rational(1, 2)`) and returns that, letting the full `Power` machinery handle all simplification (exact perfect squares, `Sqrt[8] -> 2 Sqrt[2]` radical extraction, numeric/MPFR evaluation, infinity algebra). `Sqrt` carries `LISTABLE | NUMERICFUNCTION | PROTECTED`.
