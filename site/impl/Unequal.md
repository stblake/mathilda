---
source: src/comparisons.c
---
`builtin_unequal` is the value-level (not structural) negation of `Equal`. For every argument pair it first checks `expr_eq`; on failure it tries `compare_numeric` to decide equality/inequality numerically. If any pair is found equal it returns `False` immediately. It returns `True` only when *every* pair is provably unequal — either decided by `compare_numeric` or, for non-comparable values, when both sides are distinct raw data (`is_raw_data`). If some pair is neither equal nor provably unequal (e.g. symbolic), it returns `NULL` so the call stays unevaluated. Fewer than two arguments returns `True`.
