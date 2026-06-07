---
source: src/comparisons.c
---
**Algorithm.** `builtin_equal` walks adjacent argument pairs. For each pair it first tests structural identity (`expr_eq`); if that fails it calls `compare_numeric`. `compare_numeric` does exact GMP comparison (`mpz_cmp`) when both sides are integer-like (so `10^30 == 10^30 + 1` is correctly False even past 2^53), exact `long double` cross-multiplied comparison when both are rational/integer, and otherwise a tolerance comparison on the doubles (relative tolerance `2^-46`) so machine reals that agree to ~14 digits compare equal. A pair compares equal → continue; a decidable non-equal pair (or two distinct "raw data" leaves, via `is_raw_data`) → return `False` immediately. If some pair is undecidable (symbolic), the whole call returns NULL (unevaluated). All-equal returns `True`. `Equal[]`/`Equal[x]` return `True`.

**Data structures.** Operates directly on the `Expr` argument array; numeric extraction goes through `get_numeric_value` (double + exact rational num/den + exactness flag) and GMP `mpz_t` for big integers.
