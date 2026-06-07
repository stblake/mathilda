---
source: src/boolean.c
---
**Algorithm.** `builtin_conditional_expression` takes `ConditionalExpression[expr, cond]` (two args; `ATTR_PROTECTED`, no Hold, so both are pre-evaluated). When `cond` is the interned `True` it yields `expr` (stealing the slot via `args[0] = NULL` so the evaluator's free of `res` doesn't double-free); when `cond` is `False` it yields the symbol `Undefined`. Nested forms `ConditionalExpression[ConditionalExpression[e, c1], c2]` are flattened to `ConditionalExpression[e, And[c1, c2]]`, with the combined `And` run through `evaluate` so contradictory or redundant conditions collapse. Any other condition leaves the call unevaluated (`NULL`).
