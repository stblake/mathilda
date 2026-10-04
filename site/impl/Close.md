---
source: src/io/streams.c
---
**Algorithm.** `builtin_close` resolves its argument two ways: an `InputStream`/`OutputStream`
object through `stream_handle` (reading the id out of the object's second slot), or a `"file"`
/`File["file"]` name through `stream_find_by_name(name, -1)`, which returns the
**most-recently-opened** matching stream. With the slot found, it captures the file name into a
fresh string, calls `stream_close_id` — which `free`s the slot's `name` and `buf`, `fclose`s
the `fp` if any, and `memset`s the slot so its `id` returns to `0` (free) — and **returns the
file name**. A stream that is not open prints `Close::stream` and returns `$Failed`.

**Data structures — the stream layer.** Operates on the process-global `Stream` registry
described under `OpenRead`. Closing is what frees a slot deterministically; the
monotonic `g_next_id` guarantees the freed id is never handed to a later open, so a stale
object cannot alias a new stream. Returning the name (rather than `Null`) is what lets `Close`
be the convenient last line of a round trip.

**Complexity / limits.** `O(n)` over the registry to find a slot by id or name (`n` = open
streams, tiny in practice). For an output stream the `fclose` is what guarantees the file is
flushed to disk. Pure ANSI C99 (`fclose`). `ATTR_PROTECTED`.
