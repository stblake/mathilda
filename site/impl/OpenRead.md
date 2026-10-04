---
source: src/io/streams.c
---
**Algorithm.** `builtin_openread` calls `open_common(res, output=0, append=0, "OpenRead")`,
which pulls the name out of a `"file"` string or `File["file"]` (`stream_filename_arg`) and
calls `stream_open_input`. That **slurps the whole file** into a `malloc`'d buffer
(`slurp_file`: `fopen("rb")`, `fseek`/`ftell` for the size, one `fread`, NUL-terminate),
claims a registry slot, and records `{name, buf, len, pos = 0}`. The handed-back value is the
inert object `InputStream["file", id]` built by `stream_make_object`, where `id` is the slot's
monotonic handle. A file that cannot be opened prints `OpenRead::noopen` and returns `$Failed`.

**Data structures — the stream layer.** `src/io/streams.c` keeps one **process-global stream
registry**: a dynamically grown array of `Stream` slots (`alloc_slot` doubles it; a slot with
`id == 0` is free), with `g_next_id` a monotonic handle allocator so a reused slot never
reuses a live id. An input stream is a resident buffer plus a moving **current point** `pos`;
the integer in `InputStream["file", id]` indexes the registry across calls, which is how
successive `Read` / `StreamPosition` calls address the same
stream. A `streams_shutdown` `atexit` hook frees every slot's `name`/`buf` and `fclose`s any
open `fp`, so nothing leaks under valgrind.

**Complexity / limits.** `O(file size)` time and memory for the initial slurp (the entire file
is resident, so a very large file is fully loaded), then `O(1)` per seek. Pure ANSI C99
(`fopen`/`fread`/`fseek`/`ftell`), so no POSIX feature-test guards are needed. `ATTR_PROTECTED`.
