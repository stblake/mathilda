---
source: src/piecewise.c
---
**Algorithm.** `builtin_round` is `do_piecewise(res, OP_ROUND, "Round", allow_2_args=true)`, sharing the rounding kernel with `Floor`/`Ceiling`/`IntegerPart`/`FractionalPart`. Round uses **banker's rounding** (round-half-to-even) everywhere. For machine `double`/`EXPR_REAL` values it calls `round_half_even`, which floors and inspects the fractional residue, breaking exact `.5` ties toward the even integer. For exact rationals it works entirely in GMP: it computes `floor((2·num + den) / (2·den))` (round-half-up) and corrects the exact-half tie (`rem == 0` and `q` odd) by decrementing `q` toward the even neighbour. MPFR values round via `mpfr_rint(..., MPFR_RNDN)`. Complex values round real and imaginary parts independently. An exact non-rational numeric quantity (e.g. `10000000 Pi`) is resolved through `do_piecewise_numeric_exact`, which numericalises at growing precision until the integer result is unambiguous.

The two-argument form `Round[x, a]` is rewritten as `a · Round[x/a]` (rounding to the nearest multiple of `a`). Symbolic arguments return `NULL`. Attributes include `ATTR_LISTABLE | ATTR_NUMERICFUNCTION`.
