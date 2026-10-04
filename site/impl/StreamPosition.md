---
source: src/io/streams.c
---
**Algorithm.** `builtin_streamposition` resolves its argument with `resolve_open` — a stream
object (`stream_handle`) or an already-open file name (`stream_find_by_name`), with **no
auto-open**: a stream that is not already open prints `StreamPosition::stream` and the call
returns `$Failed`. The position is read as a byte offset: for an output stream it is
`ftell(s->fp)`, for an input stream it is the buffer cursor `(long)s->pos`. A negative `ftell`
is clamped to `0`, and the offset is returned as an `Integer`.

**Data structures — the stream layer.** Reads the current point of a slot in the process-global
`Stream` registry described under `OpenRead`: an input stream's `pos` into its
resident buffer, or an output stream's file offset. It is the query half of the seek pair whose
setter is `SetStreamPosition`; reading through `Read`
advances this same `pos`.

**Complexity / limits.** `O(n)` over the registry to locate the slot, then `O(1)`. The offset
is in **bytes**, not characters or objects. Pure ANSI C99 (`ftell`). `ATTR_PROTECTED`.
