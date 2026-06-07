---
source: src/list.c
---
**Algorithm.** `builtin_max` flattens any `List` arguments into a flat argument sequence, then
scans for the maximum among real-numeric atoms (compared with `expr_compare`) while collecting
distinct non-numeric/symbolic terms. `Infinity`/`-Infinity` and `Overflow[]` are handled as
absorbing/identity elements. If everything reduces to numbers it returns the single largest
value; otherwise it returns `Max[...]` over the numeric maximum plus the remaining symbolic
terms (returning `NULL` to stay unevaluated when nothing simplified). Empty `Max[]` is
`-Infinity`.
