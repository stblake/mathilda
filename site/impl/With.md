---
source: src/modular.c
---
**Algorithm.** `builtin_with` (in `src/modular.c`) implements lexical constants by **literal substitution** — no renaming and no symbol-table mutation (contrast `Module`/`Block`). `With` carries `HoldAll | Protected` (set in `src/attr.c`). Each binding must be `x = val` (value evaluated immediately in the outer scope) or `x := val` (RHS substituted verbatim, unevaluated); the handler collects these into a `ScopingEnv` linked list mapping the constant name to its replacement expression.

The body is rewritten by `substitute_scoping`, the same recursive, shadow-aware tree walk used by `Module`: it replaces free occurrences of each constant, drops a name from the environment when descending into a nested scoping construct that rebinds it, and substitutes into nested binding RHSs (so `With[{q=12}, With[{k=q}, k]]` resolves `k` to 12) without touching binding LHS names. The substituted body is then `evaluate`d; `Return[v]`/`Return[v, With]` is trapped via `eval_classify_return`. Because substitution is structural, the constants vanish before evaluation and leave no trace in the symbol table.
