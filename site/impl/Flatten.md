---
source: src/list.c
---
**Algorithm.** `builtin_flatten` (in `src/list.c`) accepts `Flatten[list]`, `Flatten[list, n]` (level cap), and `Flatten[list, n, h]` (custom head). It iterates over the top-level arguments calling the recursive worker `flatten_rec`, which splices the children of any subexpression whose head equals the flattening head `h` (default `List`) up into the output, descending up to `n` levels (n = -1 means unlimited). The collected arguments are gathered into a growable buffer and reassembled under the original head.

**Data structures.** A dynamically grown `Expr**` accumulator (`results`, with `count`/`cap`) holds the deep-copied leaf expressions before the final `expr_new_function` rebuild.
