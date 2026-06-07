---
source: src/solve.c
---
`InverseFunctions` is an option symbol for `Solve`, not a callable function. It is interned in `sym_names.c` and documented in `solve.c`; it has no builtin handler. `solve.c`'s option parser recognises `InverseFunctions -> spec`: the default `Automatic` (and any value other than `False`) leaves the `enabled` flag set in the solver state (`solveinv.h`), while `InverseFunctions -> False` disables it. When enabled, `Solve` is allowed to "peel" invertible heads (Exp, Log, trig, powers) off an equation by applying their inverse function — the work is done by the inverse-function specialist registered as `` Solve`SolveInverseFunctions `` in `solveinv.c`. The symbol is consumed only as a configuration key.
