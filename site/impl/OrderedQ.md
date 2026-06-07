---
source: src/sort.c
---
**Algorithm.** `builtin_orderedq` scans adjacent element pairs of the list once and returns
`True` iff every pair is in canonical order. With no ordering function it uses `expr_compare`
(the canonical order described under `Sort`); with a second argument `p` it evaluates `p[a, b]`
and treats `True`/`1` as ordered. The first out-of-order pair short-circuits to `False`. Empty
and single-element lists are trivially `True`. This is the linear-scan predicate corresponding
to the same comparator that `Sort` uses.
