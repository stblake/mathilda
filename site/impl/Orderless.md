---
source: src/attr.c
---
`Orderless` is an attribute marker symbol, not a function. `attr.c` maps the name to the `ATTR_ORDERLESS` bitflag. A head carrying this bit makes the evaluator's ordering step (`eval_sort_args` in `eval.c`) sort the arguments into canonical order (`expr_compare`), giving commutative behaviour. The pattern matcher also accounts for `ATTR_ORDERLESS` so a pattern can match commuted arguments. Plus and Times set this bit. The symbol appears only as an argument to `Attributes`/`SetAttributes`, or inside a `Function[..., Orderless]` attribute spec (`purefunc.c` maps `SYM_Orderless` → `ATTR_ORDERLESS`).
