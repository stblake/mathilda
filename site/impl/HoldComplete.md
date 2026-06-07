---
source: src/attr.c
---
`HoldComplete` has no C handler; the attribute table in `src/attr.c` gives it `ATTR_HOLDALLCOMPLETE | ATTR_PROTECTED`, so the evaluator suppresses all argument evaluation and the upvalue/Sequence/Unevaluated-stripping machinery as well. `ReleaseHold` removes the wrapper.
