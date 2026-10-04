---
source: src/message.c
---
**Algorithm.** `Message` is `HoldFirst, Protected`: the message name `sym::tag`
(a `MessageName`) is held, so that evaluating it yields the template string the
user defined for it, while the trailing arguments evaluate normally.
`builtin_message` first calls `mth_msg_note_fired()` unconditionally — so an
enclosing `Check` registers the diagnostic whether or not it is displayed — then,
unless messages are suppressed, evaluates the first argument and, when it resolves
to a string, prints it to stderr. It always returns `Null`.

**Data structures.** It shares the message subsystem's single fired-counter and
suppression depth with `Check` and `Quiet`; there is no per-message table beyond
the `MessageName` own-value that holds a template string.

**Complexity / limits.** Argument substitution into the template (the `` `1` ``
slots) is **not** performed: the port uses `Message` only in error branches that
are immediately followed by `Throw`, so the printed text is diagnostic rather than
load-bearing. The firing is what matters — it is the hook that `Check` detects and
`Quiet` silences, both through the same funnel.
