---
source: src/iter.c
---
**Algorithm.** `Break[]` is a zero-argument flow-control marker, `Protected` with
no `Hold` attributes. `builtin_break` only validates arity — a non-zero argument
count emits the standard `Break::argx` message through `builtin_arg_error` and
leaves the call unevaluated — and otherwise returns `NULL`, so the raw `Break[]`
node stands unevaluated as the marker itself. It belongs to the **head-detected**
family of non-local control (Mechanism B: `Return`/`Break`/`Continue`/`Abort`/
`Quit`), distinct from the `Throw`/`Goto` sentinel. The difference is where it is
seen: `Do`/`For`/`While` (`src/iter.c`) inspect it by head through
`iter_flow_classify` (keyed on the interned `SYM_Break`) at the loop boundary,
but it is **not** short-circuited in `evaluate_step`'s argument-evaluation loop,
so `Print[Break[]]` does not escape — the semantics are lexical.

**Data structures.** None; `Break[]` carries no payload and is just an unevaluated
`EXPR_FUNCTION` node that the loops match by pointer-identity of its head symbol.

**Complexity / limits.** A `Break[]` takes effect as soon as it is evaluated and
escapes only the *innermost* enclosing `Do`/`For`/`While`, which then yields
`Null`. `Table` deliberately does not honour it. A `Break[]` that reaches top
level with no loop to consume it is reported by
`eval_report_uncaught_break_continue` (`src/eval.c`) with `Break::nofwd` and
rewritten to the inert `Hold[Break[]]`, so feeding it back does not re-trigger.
