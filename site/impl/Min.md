---
source: src/list.c
---
**Algorithm.** `builtin_min` mirrors `Max`: it flattens `List` arguments, scans real-numeric
atoms for the minimum (via `expr_compare`), collects distinct symbolic terms, and treats
`Infinity`/`-Infinity`/`Overflow[]` as identity/absorbing elements. All-numeric input returns
the single smallest value; mixed input returns `Min[...]` over the numeric minimum and the
remaining symbolic terms, or `NULL` if nothing simplified. Empty `Min[]` is `Infinity`.
