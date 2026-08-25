---
source: src/core.c
---
`builtin_quotientremainder` returns the pair `{Quotient[m, n], Mod[m, n]}`, sharing the floored-division and Gaussian-integer conventions of its two components, so the remainder always carries the sign of the divisor. Registered `PROTECTED | NUMERICFUNCTION | LISTABLE`; non-numeric arguments are left unevaluated.
