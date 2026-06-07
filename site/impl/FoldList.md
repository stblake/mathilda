---
source: src/funcprog.c
---
`builtin_foldlist` is `fold_impl(res, true)`: it returns every intermediate
accumulator, `FoldList[f, x, {a,b,c}]` → `{x, f[x,a], f[f[x,a],b],
f[f[f[x,a],b],c]}`. The two-argument form seeds with the list's first element;
`FoldList[f, {}]` returns an empty list under the input list's head. The seed is
pushed into an `ExprBuf`, and `iter_run` with `fold_step` consumes the remaining
elements left-to-right via `apply_binary(f, acc, elem)` (`eval_and_free`). With
`as_list=true`, `ebuf_finalize` wraps the full history under the **input list's
head** (preserved via `expr_copy(list_head)`). Shares all machinery with `Fold`.
