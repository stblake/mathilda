---
source: src/list.c
---
`builtin_total` sums the elements of a list, optionally restricted to a level range. `Total[list]` sums the top level; `Total[list, n]` sums levels 1..n; `Total[list, {n}]` sums exactly level n; `Total[list, {n1, n2}]` sums a range; `Total[list, Infinity]` sums all levels. Negative level indices count from the bottom using the list's depth (`get_depth_for_total`). The chosen levels' elements are gathered and combined with `Plus` (so the usual numeric/symbolic Plus folding applies — `Total` is just structural element collection feeding `Plus`). `Total` carries `ATTR_PROTECTED`.
