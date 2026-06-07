---
source: src/list.c
---
**Algorithm.** `builtin_deleteduplicates` (in `src/list.c`) keeps the first occurrence of each distinct element. The default (no test) path builds a `HashTable` keyed by expression hash/equality and inserts each element only if `ht_find` reports it absent, giving expected O(N) behavior. When a custom equivalence test is supplied as the second argument it falls back to an O(N²) scan, evaluating `test[elem, kept_j]` and treating a `True` result as a duplicate.

**Data structures.** The `HashTable` from `src/list.c`'s hashing utilities (open-addressing over `Expr*` via `expr_hash`/`expr_eq`); a pre-sized `Expr**` collects the survivors before the result `List` (the original head) is built.
