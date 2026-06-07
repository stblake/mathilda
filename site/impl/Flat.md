---
source: src/attr.c
---
`Flat` is an attribute marker symbol, not a function. `attr.c` maps the name to the `ATTR_FLAT` bitflag (via `SetAttributes`/the attribute name table). When a head carries `ATTR_FLAT`, the evaluator's flattening step (`eval_flatten_args` in `eval.c`) splices nested same-head calls into the parent argument list, giving associative behaviour (`f[a, f[b, c]]` → `f[a, b, c]`). The pattern matcher (`match.c`) also consults `ATTR_FLAT` when matching sequence patterns against such heads. Plus and Times set this bit. The symbol `Flat` itself only ever appears as an argument to `Attributes`/`SetAttributes` or inside a `Function[..., Flat]` attribute spec (`pure_function_attributes` in `purefunc.c` maps `SYM_Flat` → `ATTR_FLAT`).
