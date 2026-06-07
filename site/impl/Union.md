---
source: src/list.c
---
**Algorithm.** `builtin_union` concatenates the elements of all argument lists (which must
share a common head), sorts the combined `Expr**` array with `qsort` under the canonical
`expr_compare` order, then removes adjacent duplicates — `expr_eq` by default, or an optional
`SameTest -> f` which is evaluated per adjacent pair. The result is the sorted, deduplicated
list. (`DeleteDuplicates` in the same file does the order-preserving variant using a hash table
keyed on `expr_hash`/`expr_eq`.)
