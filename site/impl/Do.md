---
source: src/iter.c
---
**Algorithm.** `builtin_do` is `ATTR_HOLDALL`, so the body and iterator spec are re-evaluated on every pass rather than once at call time. A multi-spec `Do[expr, s1,...,sk]` is rewritten as nested two-spec `Do[Do[expr, sk], s1,...,s_{k-1}]` and handed back to the evaluator, reducing every case to a single iterator. The spec is parsed by the shared `iter_spec_parse` into an `IterSpec` (kinds `COUNT` for `n`/`{n}`, `RANGE` for `{i,imin,imax,di}` with defaults, `LIST` for `{i,{...}}`); numeric bounds are resolved to doubles by `iter_spec_resolve_numeric`, with `Infinity` allowed for unbounded loops.

**Variable localization.** Before looping, `iter_spec_shadow` saves and clears the iterator symbol's `own_values`; each iteration binds the current value via `symtab_add_own_value`, and `iter_spec_restore` frees the per-iteration binding chain and restores the original OwnValue afterward — a manual mimic of Mathematica's iterator scoping. For exact ranges the bound value `curr_e` is advanced with `Plus[curr_e, di_e]` evaluated each step (keeping integer/rational exactness), while a parallel `double val` drives the termination comparison (with a `1e-14` tolerance).

**Control flow.** After each `evaluate(body)`, `iter_flow_classify` (boundary head `SYM_Do`) maps the result to break / continue / return-value / propagate (Throw/Abort/Quit/foreign Return). `Continue` in a range loop still advances the counter before re-testing. Returns the Return payload if any, else `Null`.
