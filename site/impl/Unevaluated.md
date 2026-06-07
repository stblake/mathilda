---
source: src/eval.c
---
`Unevaluated` has no C handler; the attribute table (`src/attr.c`) gives it `ATTR_HOLDALLCOMPLETE | ATTR_PROTECTED`. The evaluator (`src/eval.c`) special-cases it: in a non-held argument position `f[Unevaluated[expr]]` passes `expr` itself (unevaluated) to `f`, stripping the wrapper. In held positions and under `HoldAllComplete` heads the wrapper is preserved.
