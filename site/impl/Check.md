---
source: src/message.c
---
**Algorithm.** `Check` is `HoldAll, Protected`, so its argument arrives
unevaluated and is evaluated under `Check`'s watch. `builtin_check` snapshots the
global message-fired counter (`mth_msg_fired_count`), evaluates `expr`, and
compares the counter afterward: if any diagnostic fired during the evaluation it
frees `expr`'s value and returns the evaluated `failexpr`, otherwise it returns
`expr`'s value. The optional third `spec` argument is accepted and ignored — any
message counts.

**Data structures.** A single `unsigned long` counter incremented by
`mth_msg_note_fired()` from inside the message funnel (`mth_message_v`). Crucially
the funnel notes a firing *even while `Quiet[]` is suppressing the print*, which
is what makes the `Quiet[Check[expr, failexpr]]` idiom — detect a failure without
showing its message — work. This is why every user-facing diagnostic must route
through the funnel: a raw `fprintf(stderr, …)` would be invisible to this counter
and `Check` would silently take the wrong branch.

**Complexity / limits.** `O(1)` bookkeeping around the inner evaluation. A `Throw`
inside `expr` is not a message: `builtin_check` tests `eval_is_inflight_throw` on
the value first and lets the sentinel propagate unchanged, so a throw escapes
`Check` rather than being reported as a failure.
