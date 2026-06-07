---
source: src/solve.c
---
`GeneratedParameters` is an option symbol for `Solve`, not a callable function. It is interned in `sym_names.c` and given a docstring in `solve.c`; it has no builtin handler. `solve.c`'s option parser recognises `GeneratedParameters -> head` and stores the (interned) head string in the solver state's `param_head` field (declared in `solveinv.h`, default `"C"`). When `Solve` returns a parameterised family (e.g. integer-multiple solutions of a trig equation), the inverse-function specialist (`solveinv.c`) emits the free parameters as `head[k]` — so `GeneratedParameters -> C` produces `C[1]`, `C[2]`, … This is purely a configuration key consulted during solving.
