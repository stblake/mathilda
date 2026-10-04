---
source: src/message.c
---
**Algorithm.** `Quiet` is `HoldAll, Protected`, so the argument is evaluated
*under* the suppression rather than before it. `builtin_quiet` brackets the
evaluation between `mth_msg_suppress_push()` and `mth_msg_suppress_pop()` and
returns the value (an in-flight `Throw` sentinel propagates unchanged). The
optional `spec` of the two-argument form — a message name or list — is accepted
and ignored: all messages are suppressed, a harmless superset of any requested
set.

**Data structures.** A single static nesting depth `g_msg_suppress_depth`
(`> 0` ⇒ suppress). The depth is a *count*, so `Quiet` nests correctly. Inside the
funnel `mth_message_v` still calls `mth_msg_note_fired()` before testing
`mth_msg_suppressed()`, so a message under `Quiet` still *fires* (an enclosing
`Check` sees it) and only its printing is silenced. A paired
`mth_msg_suppress_depth_save`/`_load` lets `TimeConstrained`'s `siglongjmp` restore
the depth if a timeout unwinds out of a `Quiet` region, which would otherwise
leave messages silenced for the rest of the session.

**Complexity / limits.** `O(1)` push/pop around the inner evaluation. `Quiet`
affects only the *display* of diagnostics; it changes neither the value computed
nor whether a message is considered to have fired.
