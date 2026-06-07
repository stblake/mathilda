---
source: src/list.c
---
**Algorithm.** `builtin_split` partitions the list into maximal runs of consecutive equal
elements, each run wrapped in the list's head. It scans left to right comparing each element to
its predecessor: `expr_eq` by default, or an optional two-argument test evaluated as
`test[prev, curr]`. A run boundary is placed wherever the test fails (and at the end), and the
accumulated run is emitted. Output is a list of run sublists.
