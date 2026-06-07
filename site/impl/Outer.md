---
source: src/funcprog.c
---
**Algorithm.** `builtin_outer` computes the generalised outer product `Outer[f, t1, t2, ..., {n1, ...}]`. It first counts trailing Integer / `Infinity` arguments as per-tensor depth limits (default `INT64_MAX`, i.e. descend to the leaves), leaving the remaining arguments after `f` as the input tensors. The recursive worker `outer_rec` walks each tensor down to its target depth, collecting one atom per tensor into `current_atoms`, and at the deepest level emits `f[a1, a2, ...]`; the assembled tree is then `evaluate`-d once. The result head for the assembled levels is taken from the first function-typed tensor.

**Data structures / limits.** Nested `Expr*` `List`s; per-tensor depths in an `int64_t[]`. With no tensors, `f[]` is returned evaluated. This is the generic functional-programming `Outer`, not a linear-algebra-specific kernel; `KroneckerProduct` and matrix outer products are built on it.
