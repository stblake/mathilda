---
source: src/core.c
---
`builtin_evaluate` (`src/core.c`) returns a copy of its single argument. Its real effect happens earlier: the evaluator forces `Evaluate[expr]` arguments to be evaluated even inside a `Hold*` head's held positions, so by the time the builtin runs the argument is already evaluated and it merely unwraps it.
