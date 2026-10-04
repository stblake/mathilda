---
source: src/funcprog.c
---
**Algorithm.** `Throw` is `Protected` with no `Hold` attributes, so `value`,
`tag` and `f` are evaluated by the ordinary argument loop before the throw
propagates. `builtin_throw` does almost nothing: it validates arity (1–3) and
returns `NULL`, because the plain `Throw[...]` node **is** the in-flight sentinel.
`eval_is_inflight_throw` recognises it by head (`SYM_Throw`, arity 1–3), and
`evaluate_step`'s argument-evaluation loop short-circuits when an evaluated
argument is such a sentinel — freeing the sibling arguments and the head and
returning the sentinel up the normal return path — so a `Throw` anywhere inside a
surrounding expression unwinds to the nearest enclosing `Catch`.

**Data structures.** No separate payload type: the sentinel is the user-visible
`Throw[value]` / `Throw[value, tag]` / `Throw[value, tag, f]` tree. Propagation is
by return value, never `longjmp`, so per-frame cleanup always runs (leak-free).

**Complexity / limits.** Unlike `Return`, which only escapes a scope boundary, a
`Throw` passes through *any* enclosing head. A tagged throw is caught only by a
`Catch[expr, form]` whose `form` matches the (re-evaluated) tag. If it reaches top
level uncaught, `eval_report_uncaught_throw` (`evaluate()`) emits `Throw::nocatch`
and returns `Hold[Throw[...]]`, except that an uncaught `Throw[value, tag, f]`
returns `f[value, tag]`.
