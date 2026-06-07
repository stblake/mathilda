---
source: src/arithmetic.c
---
`Rational[n, d]` is the internal head for exact rationals. `builtin_rational` only fires for two integer arguments: it calls `make_rational(n, d)` to reduce to lowest terms with a positive denominator. If the input is already in canonical form (no reduction happened) it returns `NULL` so the structural `Rational[n, d]` is left as-is; otherwise it returns the reduced form (an `EXPR_INTEGER` when the denominator becomes 1). Division by zero emits `Power::infy` and returns `ComplexInfinity` (or `Indeterminate` for `0/0`). Non-integer arguments return `NULL`.
