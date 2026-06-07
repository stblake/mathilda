---
source: src/linalg/tr.c
---
**Algorithm.** `builtin_tr` is the generalised trace `Tr[list, f, n]`. It walks the diagonal of a rank-`n` nested `List`: `extract_diagonal_element` descends `n` levels always taking element `index` at each level, yielding the `i`-th diagonal entry, and the loop collects entries until an index falls out of bounds. The collected leaves are combined with the head `f` (default `Plus`) via a single `f[d_0, d_1, ...]` call run through the evaluator (`eval_and_free`), so for a 2-D matrix this is the ordinary sum of diagonal entries.

**Data structures / limits.** A growable `Expr**` of the diagonal leaves. The depth `n` defaults to the nesting depth of the leftmost spine (`get_default_trace_depth`, min 1) and may be given explicitly as a non-negative Integer. A non-`List` argument is returned unchanged; a non-integer explicit depth leaves the call unevaluated.
