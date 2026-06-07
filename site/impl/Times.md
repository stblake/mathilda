---
source: src/times.c
---
**Algorithm.** `Times` carries `FLAT | ORDERLESS | LISTABLE | NUMERICFUNCTION | ONEIDENTITY`, so the evaluator flattens nested `Times`, threads over `List` args, and sorts canonically before `builtin_times` runs. The builtin does numeric folding and exponent collection.

It accumulates the rational/numeric coefficient into a running `num_prod` (with inexact contagion when a Real/MPFR factor is present); a numeric factor that cannot fold (e.g. an integer carrying a pending radical) is stashed as its own group. Each remaining factor is split into a `(base, exponent)` pair (a `BasePower`; bare factors get exponent 1, `Power[b,e]` is unpacked), and factors with structurally equal base are merged by summing exponents (`Plus` over the two exponents). A radical-canonicalisation pass folds factors of a positive integer base `b ≥ 2` appearing in the accumulated rational coefficient into a `Power[b, q]` group with rational q. The result is rebuilt as the numeric coefficient times the surviving `Power` groups; `0` short-circuits and a zero result returns 0.

**Data structures.** The `BasePower` struct (`{base, exponent}` — shared with `trig_canon.h`) in a linear `groups[]` array; coefficient arithmetic on exact int/rational/bigint with `eval_and_free` for symbolic recombination. Like-base lookup is linear scan with `expr_eq`. An overflow guard returns an `Overflow[]` sentinel.

**Complexity / limits.** Roughly `O(n^2)` from the linear base-grouping over n factors (canonical sort tends to make equal bases adjacent).
