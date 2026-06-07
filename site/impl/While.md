---
source: src/iter.c
---
**Algorithm.** `builtin_while` implements `While[test]` or `While[test, body]` and is `ATTR_HOLDALL`, so both arguments are re-evaluated each iteration. Each pass calls `evaluate(test)`; the result is first run through `iter_flow_classify` (boundary head `SYM_While`) so flow-control raised inside the test itself is honored (Break/Return exit, Throw/Abort/Quit propagate, Continue restarts the loop). The loop then exits unless the test is exactly the symbol `True`. The body (if present) is evaluated and classified the same way: `Break` exits (yielding `Null`), `Continue` skips to the next test, `Return[v]` exits yielding `v`. Truth testing is pointer equality on `SYM_True`. The loop returns the Return payload if one was issued, else `Null`.
