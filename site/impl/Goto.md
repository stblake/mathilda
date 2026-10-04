---
source: src/funcprog.c
---
**Algorithm.** `Goto` is `Protected` with no `Hold` attributes, so `tag` is
evaluated by the argument loop. `builtin_goto` validates arity (exactly one
argument) and returns `NULL`, because — exactly as with `Throw` — the plain
`Goto[tag]` node is itself the in-flight sentinel. `eval_is_inflight_goto`
recognises it by head, and `evaluate_step`'s argument loop short-circuits on it,
so a `Goto` fired inside a nested call (an `If` branch, say) bubbles up to the
enclosing `CompoundExpression`. That `CompoundExpression` (`builtin_compound`
expression) is the construct that actually consumes the sentinel: it scans its own
statements for a raw held `Label[tag]` whose tag is structurally equal and resumes
evaluation there — a forward jump skipping intervening statements, or a backward
jump forming a loop. If the current compound expression has no matching `Label`,
the sentinel propagates to the enclosing one.

**Data structures.** The sentinel is the `Goto[tag]` node itself; targets are raw
held `Label[tag]` nodes scanned linearly. Propagation is by the normal return path
(no `setjmp`/`longjmp`), so it is leak-free.

**Complexity / limits.** A `Goto` loop is a genuine loop with no artificial
iteration cap; termination is the program's responsibility, as with `While`. A
`Goto[tag]` that reaches top level with no matching `Label` anywhere emits
`Goto::nolabel` (stderr) and returns the inert `Goto[tag]` node; the message fires
only when truly unmatched, so a `Goto` that legitimately propagates from an inner
to an outer `CompoundExpression` mid-evaluation stays silent.
