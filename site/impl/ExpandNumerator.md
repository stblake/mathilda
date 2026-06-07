---
source: src/expand.c
---
**Algorithm.** `builtin_expand_numerator` (in `src/expand.c`) calls `expr_expand_numerator`, which separates an expression's numerator from its denominator and expands only the numerator. For a `Times`, it partitions factors into denominator factors (those of the form `Power[base, k]` with `k` a negative integer, detected by `is_negative_int_power`) and the rest; the non-denominator product is run through `expr_expand`, then recombined with the untouched denominator factors. A bare `Power` with negative integer exponent is a pure denominator and is returned unchanged; a positive/symbolic power is expanded at the top level. It threads over `List`, equations, inequalities, `And`/`Or`/`Not`, and `Plus` (the `is_thread_head` set), expanding per-summand.

**Data structures.** Separate `Expr**` accumulators for numerator and denominator factors; the denominator is preserved verbatim while only the numerator passes through `expr_expand`.
