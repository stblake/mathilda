---
source: src/trace.c
---
**Algorithm.** `builtin_trace` (`src/trace.c`) is a thin wrapper over the
evaluator-side collector `eval_collect_trace` (`src/eval.c`), which re-runs the
argument's evaluation while recording every form it passes through. Because
`Trace` carries `HoldAll`, the argument reaches the builtin unevaluated, so its
rewrite sequence is observed from the start (and the evaluation clock is bumped
once so an already-evaluated argument is still traced in full). The collector
returns a **nested** `List` mirroring the evaluator's own recursion: each argument
or head sub-evaluation that takes a step becomes a sublist, one that takes no step
contributes nothing, and the reassembled intermediate form appears as a step.
`Trace[expr, form]` then filters that tree with `trace_collect_matches`, keeping
only step leaves that structurally `match` the held pattern `form` and flattening
the nesting into a plain `List`.

**Data structures.** The raw collector output is an `Expr` tree whose `List` nodes
are the nesting markers. `trace_holdform_tree` walks it and wraps every non-`List`
leaf in `HoldForm`, so the returned structure is inert under the evaluator's
fixed-point re-pass (a bare `1 + 1` would otherwise reduce again to give `{2, 2}`)
while still printing transparently. The two-argument filter accumulates matches in
a growable `Expr**` buffer, copying each matched node.

**Complexity / limits.** Proportional to the number of evaluation steps `expr`
takes. A builtin's internal computation and `Listable` threading show as a single
atomic rewrite (matching Mathematica — `Range[10]` is one step). `Trace` is
reentrant (an inner `Trace` appears as one reduced value to the outer). Arities
other than 1 or 2 return `NULL` (stay unevaluated); `TraceDepth` is not
implemented. Attributes `HoldAll`, `Protected`.
