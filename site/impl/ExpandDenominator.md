---
source: src/expand.c
---
**Algorithm.** `builtin_expand_denominator` (in `src/expand.c`) calls `expr_expand_denominator`, the mirror of `ExpandNumerator`: it expands only the denominator. For a `Times` it collects the negative-integer-power factors (`is_negative_int_power`), rebuilds them as positive powers, multiplies them into a single denominator product, runs `expr_expand` on that product, and re-inverts it as `Power[expandedDen, -1]` multiplied back against the untouched numerator factors. A bare `Power[base, -k]` is handled by expanding `base^k` and re-inverting. Threading over `List`, equations, inequalities, logic heads and `Plus` matches `ExpandNumerator`.

**Data structures.** Parallel `Expr**` buffers for numerator factors and positive-power denominator factors; only the denominator product is expanded, the numerator is copied through.
