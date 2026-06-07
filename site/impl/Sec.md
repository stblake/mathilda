---
source: src/trig.c
---
**Algorithm.** `builtin_sec` follows the `src/trig.c` cascade but uses even symmetry: `strip_inverse_call(arg, "ArcSec")` for `Sec[ArcSec[x]] -> x`; `even_fold` for `Sec[-x] -> Sec[x]` when the argument is superficially negative; `trig_i_fold(arg, "Sech", 0)` for `Sec[I y] -> Sech[y]`; and `Sec[0] = 1`. Exact values at rational multiples of `Pi` are recognised by `extract_pi_multiplier` and produced by `exact_sec`.

**Numeric.** MPFR arguments use `numeric_mpfr_apply_unary(..., mpfr_sec)` (complex fallback `mpfr_complex_sec`); otherwise `get_approx` computes `1.0 / ccos(c)` for inexact inputs, yielding `EXPR_REAL` or `Complex`. Symbolic input returns `NULL`. Attributes: `ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED`.
