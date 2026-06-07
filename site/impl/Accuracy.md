---
source: src/precision.c
---
`builtin_accuracy` (`src/precision.c`) delegates to `accuracy_of`, which returns the number of correct digits to the right of the decimal point. Exact quantities (integers, bigints, exact rationals, strings, symbols, exact zero) return `Infinity`. For inexact reals it returns `MachinePrecisionDigits - log10|x|`; for `EXPR_MPFR` it uses the value's actual precision (`mpfr_get_prec / log2(10)`) minus `log10|x|` (inexact zero gets the precision in digits directly). `Complex[re, im]` and general function arguments recurse and take the minimum (`precision_min`) over components/arguments. `Accuracy` is `ATTR_LISTABLE`, so it threads over lists.
