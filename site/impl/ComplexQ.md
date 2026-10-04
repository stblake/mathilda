---
source: src/core.c
---
**Algorithm.** `builtin_complexq` is a one-argument syntactic predicate. It tests
whether its (already-evaluated) argument has head `Complex` — the single check
`core_head_is(arg, SYM_Complex)` on an interned-pointer comparison — and returns
the symbol `True` or `False` accordingly. A purely real number is *not* `Complex`:
Mathilda stores `2 + 3 I` as `Complex[2, 3]` but keeps a real as a bare `Real`, so
`ComplexQ[2.0]` is `False`. The test is on the head alone, so it does not care
whether the parts are exact or inexact.

**Data structures.** None beyond the argument node. The comparison is a pointer
equality against the interned `SYM_Complex` name; the result is a fresh
`EXPR_SYMBOL` (`True`/`False`).

**Complexity / limits.** `O(1)`. With fewer or more than one argument it declines
(returns `NULL`) and stays symbolic. Being a syntactic head test, it reports on
representation, not mathematical identity — a value that happens to be real-valued
but is written with a `Complex` head is a matter for the arithmetic that built it,
not for `ComplexQ`.
