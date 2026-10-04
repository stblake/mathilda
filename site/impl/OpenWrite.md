---
source: src/io/streams.c
---
**Algorithm.** `builtin_openwrite` calls `open_common(res, output=1, append=0, "OpenWrite")`,
which resolves the name (`stream_filename_arg` accepts `"file"` or `File["file"]`) and calls
`stream_open_output(name, append=0)` — `fopen(name, "wb")`, which **truncates** any existing
file. A fresh registry slot records `{name, fp, is_output = 1}` and the result is the inert
object `OutputStream["file", id]`. A file that cannot be opened prints `OpenWrite::noopen` and
returns `$Failed`.

**Data structures — the stream layer.** Shares the process-global `Stream` registry described
under `OpenRead`. An output stream holds an open `FILE*` rather than a buffer;
`id` in `OutputStream["file", id]` indexes the slot so `Write` /
`StreamPosition` address it across calls. The `atexit` hook `fclose`s any
stream still open at exit, so an output file is flushed even if `Close` is
forgotten — though `Write`/`WriteString` `fflush` after every
call, so a round trip with `Read` sees the bytes immediately.

**Complexity / limits.** `O(1)`; the open is a single `fopen`. The truncation is the
distinction from `OpenAppend` (`"ab"`). Pure ANSI C99 (`fopen`), no POSIX
guards. `ATTR_PROTECTED`.
