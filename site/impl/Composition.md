---
source: src/core.c
---
`builtin_composition` (`src/core.c`) handles only the algebraic simplifications of `Composition[f1,...,fn]`: `Composition[]` -> `Identity`, `Composition[f]` -> `f`, dropping `Identity` arguments, and cancelling adjacent `f`/`InverseFunction[f]` pairs (in either order, iterated to fixed point). The actual application `Composition[f1,...,fn][args]` -> `f1[f2[...fn[args]...]]` is performed in the evaluator (`eval.c`). The symbol carries `ATTR_FLAT | ATTR_ONEIDENTITY`.
