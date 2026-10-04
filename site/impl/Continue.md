---
source: src/iter.c
---
**Algorithm.** `Continue[]` is a zero-argument flow-control marker, `Protected`
with no `Hold` attributes. `builtin_continue` validates arity only — a non-zero
argument count emits `Continue::argx` through `builtin_arg_error` — and otherwise
returns `NULL`, so the raw `Continue[]` node stands as the marker. Like `Break`,
it is a **head-detected** marker (Mechanism B), recognised by the loop builtins
through `iter_flow_classify` (interned `SYM_Continue`) and not short-circuited in
the argument-evaluation loop, so it only takes effect at a loop boundary.

**Data structures.** None; a bare unevaluated node matched by head.

**Complexity / limits.** `Continue[]` skips the remainder of the current body of
the innermost `Do`/`For`/`While` and advances to the next iteration: in `Do` it
advances the iterator and re-tests (the arithmetic-progression form must still
step its running value on `ITER_FLOW_CONTINUE`, else the loop would spin), in
`For` it evaluates the increment step and re-tests, and in `While` it
re-evaluates the test. A `Continue[]` that reaches top level is reported with
`Continue::nofwd` by `eval_report_uncaught_break_continue` and rewritten to the
inert `Hold[Continue[]]`.
