---
source: src/io/streams.c
---
**Algorithm.** `builtin_streams` walks the process-global stream registry and builds a `List`
of the live stream objects. With one argument it reads a file-name filter
(`stream_filename_arg`) and keeps only slots whose `name` matches it; with none it keeps every
open slot. For each surviving slot it emits `stream_make_object` (`InputStream["file", id]` or
`OutputStream["file", id]` by the slot's `is_output` flag) into a `malloc`'d array that is
grown as needed, then wraps it in `List[...]`.

**Data structures — the stream layer.** Reads the same process-global `Stream` registry
described under `OpenRead`; the result is a snapshot of the slots whose
`id != 0`, so it shrinks as `Close` frees slots and is empty at program start and
once everything is closed. Each entry is the inert object itself, so it can be passed straight
back to `Read`/`Write`/`Close`.

**Complexity / limits.** `O(n)` over the registry, `n` the number of open slots. The id inside
each object is an internal handle whose value depends on how many streams have been opened, so
code should test membership rather than match a printed form. `ATTR_PROTECTED`.
