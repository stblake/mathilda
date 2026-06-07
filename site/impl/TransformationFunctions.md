---
source: src/simp/simp_builtins.c
---
**Algorithm.** `TransformationFunctions` is an option symbol for Simplify (it has
no builtin handler of its own). `builtin_simplify` detects
`Rule[TransformationFunctions, spec]` among its arguments and resolves it into a
`(use_builtin, user_funcs[])` pair: `Automatic` (the default) keeps the built-in
transform pipeline only; `{f1, ...}` suppresses the built-ins and uses only the
`fi`; `{Automatic, f1, ...}` runs the built-in pipeline *and* the `fi`; a bare
`f` is treated as the single-function list `{f}`. The `fi` are borrowed pointers
into the option expression (the evaluator keeps it alive across the call). After
the built-in search produces `best` (or, when built-ins are suppressed, `best =
expr_copy(input)`), `simp_apply_transformations` applies each user function to
the current best and keeps the lowest-complexity result by `score_with_func`.

**Data structures.** No state beyond the parsed option; `user_funcs` is a small
heap array of borrowed `Expr*` head expressions, freed after the search.
