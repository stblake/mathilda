---
source: src/plus.c
---
**Algorithm.** `Plus` carries `FLAT | ORDERLESS | LISTABLE | NUMERICFUNCTION | ONEIDENTITY`, so by the time `builtin_plus` runs the evaluator has already flattened nested `Plus`, threaded over `List` args, and canonically sorted the arguments. The builtin handles the numeric folding and like-term collection. `Plus[]` is 0, `Plus[x]` is x.

It first distributes `Times[-1, Plus[...]]` over the outer sum (`is_neg_of_plus`) so cancellations like `a + b - (a+b) -> 0` actually collapse. Inexact contagion (`numeric_contagion_args`) numericalises exact numeric parts in-place when any summand is an inexact Real/MPFR. The main pass splits each term into a numeric coefficient and a base via `get_coeff_base` (`c·base`, with bare numbers having a null base), then groups by structurally equal base (`TermGroup` array), summing coefficients within each group. Pure numeric terms fold into a single running total. The result is reassembled as `Plus` of `coeff·base` terms, dropping zero-coefficient groups.

**Data structures.** A linear `TermGroup[]` array of `(base, coeff)` pairs; coefficient arithmetic goes through the exact integer/rational/bigint paths and `eval_and_free` for symbolic recombination. Like-term lookup is linear scan with `expr_eq`.

**Complexity / limits.** Roughly `O(n^2)` worst case from the linear like-term grouping over n summands (after the evaluator's canonical sort, equal bases are typically already adjacent).
