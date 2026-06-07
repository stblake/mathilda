---
source: src/funcprog.c
---
**Algorithm.** `builtin_fold` (`fold_impl(res, false)`) is a left fold:
`Fold[f, x, {a,b,c}]` → `f[f[f[x,a],b],c]`. The two-argument form `Fold[f,
{a,b,c}]` uses `a` as the seed and folds over the rest; on an empty list this
form stays unevaluated. The seed is pushed into an `ExprBuf` history and the
shared `iter_run` driver is invoked with `fold_step`, which consumes the list
elements in order, computing `apply_binary(f, accumulator, elems[idx++])` (build
`f[acc, e]`, `eval_and_free`). `ebuf_finalize(..., as_list=false)` returns the
last accumulator. `Fold`/`FoldList` are the same `fold_impl`.

**Data structures.** The list's underlying argument array is borrowed (no copy of
the spine); each accumulator value is an owned `Expr*` stored in the history
buffer. The output preserves the input list's head only in the `FoldList` form.
