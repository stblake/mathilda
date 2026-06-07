---
source: src/eval.c
---
A user-visible system variable seeded in `recursion_limit_init` (`src/eval.c`) with an integer OwnValue (the default C-level limit). The evaluator enforces it with a static depth counter `g_eval_depth` incremented per `evaluate_step`; when depth exceeds the limit it aborts the call and emits `$RecursionLimit::reclim`. Assignment is special-cased in the `Set` handler: `$RecursionLimit = expr` routes through `recursion_limit_set` (called from `eval.c`'s assignment path), which rejects values below a hard floor with `$RecursionLimit::limset` and otherwise reinstalls the new integer as the symbol's OwnValue so the C counter and the symbol stay in sync.
