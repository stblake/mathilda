---
source: src/core.c
---
`builtin_numericq` (1-arg) returns `True`/`False` by calling the recursive predicate `is_numeric_quantity`. That predicate returns true for `EXPR_INTEGER`/`EXPR_REAL`/`EXPR_BIGINT`/`EXPR_MPFR`; for the named numeric constants `Pi`, `E`, `I`, `Infinity`, `ComplexInfinity`, `EulerGamma`, `GoldenRatio`, `Catalan`, `Degree`; for `Complex[...]` and `Rational[...]` heads; and for any function whose head carries `ATTR_NUMERICFUNCTION` *provided every argument is itself numeric* (recursive check). Everything else — bare symbols, non-numeric heads — yields `False`. Unlike `NumberQ`, this resolves the "would evaluate to a number" question structurally via the attribute system rather than by numericalizing.
