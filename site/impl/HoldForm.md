---
source: src/core.c
---
`HoldForm` has no C handler; it is purely an evaluation/display marker. It is given `ATTR_HOLDALL | ATTR_PROTECTED` in `core_init` (`src/core.c`) so its argument stays unevaluated, and the printer renders `HoldForm[expr]` as just `expr` (the wrapper is invisible). `ReleaseHold` strips it.
