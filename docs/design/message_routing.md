# Message routing: how a diagnostic reaches `Quiet[]` and `Check[]`

## The contract

`Quiet[expr]` suppresses the diagnostics an evaluation prints; `Check[expr,
fail]` returns `fail` when a diagnostic *fired* during `expr`. Mathilda has no
single message object and no central dispatcher — the two behaviours are a
**two-step convention** every emission site must uphold, defined in
`src/message.{c,h}`:

1. **`mth_msg_note_fired()`** — bump the global fired-counter so an enclosing
   `Check[]` (which snapshots `mth_msg_fired_count()` around its first argument)
   sees that a message happened, *even while suppressed*.
2. **guard the print** with `mth_msg_suppressed()` — the `Quiet[]` depth counter;
   while it is non-zero the text is not written.

A diagnostic emitted with a raw `fprintf(stderr, "Head::tag: …")` upholds
**neither**: `Quiet[Head[bad]]` still leaks the text, and
`Check[Head[bad], fallback]` returns `Head[bad]` instead of `fallback` — a wrong
answer, not a cosmetic slip. Before this was enforced, ~294 sites bypassed it
(e.g. `MatrixPower` escaped both; `Power::infy` was invisible to `Check`).

## The funnel

Route every user-facing diagnostic through **`mth_message(head, tag, fmt, …)`**
(in `src/message.c`). It performs both steps once:

```c
mth_message("Power", "infy", "Infinite expression 1/0 encountered.");
```

`head` and `tag` are separate arguments, so no `"::"`-bearing literal ever
reaches an `fprintf(stderr, …)` — which is exactly what lets `make
check-messages` tell a routed site from a bypassing one (see below).

Variants:

- **`mth_message_gated(extra_mute, head, tag, fmt, …)`** — for a subsystem with
  an *internal-probe* mute (`g_arith_warnings_muted` for `Power`/`Plus`/`Times`,
  `g_fm_quiet` for `FindMinimum`). Pass the **raw** flag, not `flag ||
  suppressed`. When it is set an internal probe is poking at divergent forms and
  the diagnostic is *sampling noise, not a user event*, so it is suppressed
  **entirely — neither printed nor noted**. This is why
  `Check[Limit[Sin[x]/x, x -> 0], bad]` returns `1`: `Limit` fires
  `Infinity::indet` on the `0·ComplexInfinity` probe internally, and `Check`
  must not catch it. Contrast `Quiet`, which is a *user* request to silence a
  *real* message — still noted, so `Quiet[Check[…]]` works.
- **`mth_message_v(extra_mute, head, tag, fmt, va_list)`** — the shared core;
  per-subsystem helpers (`fs_msg`, `dt_msg`, `ops_msg`, `root_warn`, `fit_warn`,
  `fm_warn`, the `N*` warn helpers, …) delegate to it instead of re-implementing
  note+guard.
- **`mth_message_cont(extra_mute, fmt, …)`** — a prefix-less continuation line
  under a diagnostic already emitted this call (no second note).

`arith_warn(msg)` in `arithmetic.h` follows the identical contract for the
fixed-string arithmetic diagnostics.

## The enforcement gate — `make check-messages`

`tools/check_message_routing.py` reads the source (comment-blanked,
multi-line-aware), finds every raw `fprintf/fputs(stderr, "Head::tag…")`
diagnostic, and diffs the set against two checked-in lists:

- **`EXEMPT`** — permanent, non-user-facing, with a reason: the funnel body
  itself (`src/message.c`); parser syntax errors (`parse.c`, before any `Check`
  frame exists); the REPL `file:line` reporter and option errors (`repl.c`); OOM
  aborts. (Debug tracing like `[voro]`/`DIUI:` and the OOM `fputs`+`abort` never
  match the `Head::tag` shape, so they need no entry.)
- **`BASELINE`** — the migration backlog, now **empty**. The gate is
  **assert-empty**: any new raw `Head::tag` stderr write fails it, and a stale
  `BASELINE` entry (a site migrated but not deleted from the list) fails it too.

Wired into the Linux CI job next to `check-c99` / `check-packed-aware`.

## Adding a new builtin

Emit diagnostics with `mth_message(head, tag, fmt, …)` (or a subsystem helper
that delegates to `mth_message_v`). Never `fprintf(stderr, "Head::tag: …")`
directly — `make check-messages` will fail. If a message genuinely cannot route
(a parser error before evaluation, an OOM abort), add it to the tool's `EXEMPT`
list with a one-line reason.
