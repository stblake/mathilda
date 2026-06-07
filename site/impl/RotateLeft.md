---
source: src/list.c
---
**Algorithm.** `builtin_rotateleft` cyclically shifts elements toward the front by `n`
(default 1) using `rotate_rec`: at each level it computes the wrapped offset `((n mod len) +
len) mod len` and reads element `i` from source index `(i + offset) mod len`. A `List`-valued
`n` applies a per-level shift amount as the recursion descends into nested lists. `RotateRight`
is implemented by negating `n` and calling the same routine.
