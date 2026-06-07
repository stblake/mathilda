---
source: src/attr.c
---
`Hold` is not a C builtin. It is registered in `attr.c`'s attribute table as `{"Hold", ATTR_HOLDALL | ATTR_PROTECTED}`, i.e. the symbol `Hold` simply carries the `ATTR_HOLDALL` attribute. When the evaluator (`eval.c`) processes `Hold[args...]` it reads that attribute before evaluating arguments and therefore leaves every argument unevaluated, returning the `Hold[...]` wrapper as-is. There is no per-head logic — holding falls entirely out of the generic attribute-driven argument-evaluation step. `HoldComplete`, `HoldPattern`, and `Unevaluated` are registered the same way with `ATTR_HOLDALLCOMPLETE`/`ATTR_HOLDALL`.
