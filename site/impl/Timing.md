---
source: src/datetime.c
---
`builtin_timing` brackets a single `evaluate(arg)` call with `clock()` (CPU time, `CLOCKS_PER_SEC`) and returns `{seconds, result}` as a two-element `List`, where `seconds` is `(end - start)/CLOCKS_PER_SEC` as an `EXPR_REAL`. It measures processor time, not wall-clock, and times a single evaluation only. (Note: the argument is evaluated explicitly inside the builtin; `Timing` is not given Hold attributes here.)
