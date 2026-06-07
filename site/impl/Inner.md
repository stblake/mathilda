---
source: src/funcprog.c
---
**Algorithm.** `builtin_inner` is the generalisation of `Dot` that replaces the elementwise multiply with an arbitrary `f` and the summation with an arbitrary `g`: `Inner[f, A, B, g]` contracts the last index of `A` with the first index of `B`, combining matched leaves with `f` and reducing each contraction list with `g` (both default chains are built so that `Inner[Times, A, B, Plus]` reproduces `Dot`). `g` defaults to `Plus` when omitted. The contraction is performed recursively by the `inner_A`/`inner_n1_A` helpers (the `n == 1` form contracts first-index-with-first-index directly), producing an unevaluated `g[f[…], …]` tree that is then run through `evaluate`. Non-`List`-structured operands, or mismatched contraction lengths, return `NULL` (leave unevaluated).

**Data structures.** Operands are walked as nested `EXPR_FUNCTION` trees keyed on `A`'s head (the result head is taken from `A`); leaves are deep-copied into `f`-applications. No flat dense buffer is used — unlike `Dot`, `Inner` works on generic heads and arbitrary combiners.
