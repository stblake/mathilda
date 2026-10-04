---
source: src/io/streams.c
---
**Algorithm.** `builtin_setstreamposition` resolves an already-open stream with `resolve_open`
(no auto-open; `$Failed` + `SetStreamPosition::stream` otherwise), then interprets its second
argument. `Infinity` means the end of the stream; otherwise the target must be an `Integer`,
and a negative value is clamped to `0`. For an **output** stream it `fseek`s (`SEEK_END` for
`Infinity`, else `SEEK_SET` to the offset) and returns `ftell`. For an **input** stream it sets
the buffer cursor `s->pos`, clamping a past-the-end target to the buffer length `s->len`
(`Infinity` → `len`). Either way it returns the resulting byte offset as an `Integer`.

**Data structures — the stream layer.** Moves the current point of a slot in the process-global
`Stream` registry described under `OpenRead`. It is the setter half of the seek
pair whose getter is `StreamPosition`; the clamping keeps the point within
`[0, len]` for an input stream, so a rewind-and-reread or skip-ahead is always in bounds.

**Complexity / limits.** `O(n)` over the registry to locate the slot, then `O(1)`. Offsets are
in **bytes**. A non-integer, non-`Infinity` target leaves the call unevaluated. Pure ANSI C99
(`fseek`/`ftell`). `ATTR_PROTECTED`.
