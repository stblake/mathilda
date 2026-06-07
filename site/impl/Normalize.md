---
source: src/linalg/normalize.c
---
**Algorithm.** `builtin_normalize` returns `expr / f[expr]`, where `f` defaults to `Norm` (which itself reduces to `Abs` for scalars). It builds `f[expr]`, evaluates it (`eval_and_free`), then returns `Times[expr, Power[norm_val, -1]]` evaluated. Because `Times` is `Listable`, the single reciprocal threads across every leaf of a vector / matrix / higher-rank tensor, and the same path handles scalars (including complex `z / Abs[z]`).

**Limits.** The zero short-circuit uses *exact* numeric-zero detection (`norm_is_numeric_zero`: literal Integer/Real/BigInt/MPFR `0`) — a zero vector is returned unchanged, but a symbolic norm that merely happens to vanish is left as a symbolic division so the input stays visible. Arity other than 1 or 2 emits `Normalize::argt`.
