---
source: src/funcprog.c
---
**Algorithm.** `builtin_thread` distributes a function over equal-length lists:
`Thread[f[{a,b},{c,d}]]` → `{f[a,c], f[b,d]}`. Arguments are `(expr, h, n)`: `h`
is the threading head (default `List`) and `n` is a position spec selecting which
of `f`'s `K` arguments take part. `thread_parse_spec` turns the spec into a
boolean mask of length `K`, handling `All`/`None`, an integer (first/last `|n|`
positions), and `{n}` / `{m,n}` / `{m,n,s}` index ranges with negative-from-end
indexing.

It then scans the masked, `h`-headed arguments to determine the common threading
length `L`; if two such arguments have different lengths the expression is
returned unchanged (a diagnostic is issued here in some systems; Mathilda elides it).
For each `k` in `0..L-1` it builds `f[...]` taking element `k` from every masked
threadable argument and copying all other arguments verbatim, then wraps the `L`
calls under `h[...]` and runs `evaluate()` so `f`'s attributes (`Listable`,
`OneIdentity`, ...) apply. Atoms and the no-threadable-argument case return a copy
of `expr`.

**Data structures.** `Expr`-tree only; a `bool* mask` of length `K`, plus a
`wrap_args` array of the `L` per-slice calls. The threading head uses the interned
`List` symbol (`expr_ref`) so `expr_eq` pointer comparisons work.
