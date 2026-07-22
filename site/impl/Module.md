---
source: src/modular.c
---
**Algorithm.** `builtin_module` (in `src/modular.c`) implements lexical scoping by **alpha-renaming** locals to unique temporaries. `Module` carries `HoldAll | Protected` (set in `src/attr.c`), so the variable list and body arrive unevaluated. The handler reads (and post-increments) the global `$ModuleNumber` counter and forms a per-invocation suffix; each local `x` (or `x = init`) becomes a fresh symbol `x$<n>` (e.g. `x$7`). Each temporary is tagged `ATTR_TEMPORARY`, and if it had an initializer (evaluated in the outer scope) that value is installed as an `OwnValue` on the renamed symbol.

The rename is applied to the body by `substitute_scoping`, a recursive tree walk over a `ScopingEnv` linked list mapping old name → replacement symbol. Crucially this walk is shadow-aware: when it descends into a *nested* scoping construct (`Module`/`Block`/`With`/`Function`/`Table`) that rebinds one of the same names, that name is dropped from the environment passed downward, so inner bindings are not corrupted; and binding RHSs are substituted with the outer environment (so `With[{q=...}, With[{k=q},...]]` resolves correctly) while binding LHS names are left intact. The renamed body is then `evaluate`d. `Return[v]` (or `Return[v, Module]`) targeting this boundary is trapped via `eval_classify_return`. Finally the `ScopingEnv`, the `VarInfo` temporaries, and the evaluated initializers are freed (the renamed symbols' OwnValues persist in the symbol table).

**Limit.** Body is taken as a single expression (arg_count must be 2).
