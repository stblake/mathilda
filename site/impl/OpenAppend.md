---
source: src/io/streams.c
---
**Algorithm.** `builtin_openappend` calls `open_common(res, output=1, append=1, "OpenAppend")`,
which resolves the name and calls `stream_open_output(name, append=1)` — `fopen(name, "ab")`.
Unlike `OpenWrite`'s `"wb"`, the `"ab"` mode **preserves** the existing
contents and positions writes at the end, and creates the file if it does not exist. A slot
records `{name, fp, is_output = 1}` and the result is `OutputStream["file", id]`. A file that
cannot be opened prints `OpenAppend::noopen` and returns `$Failed`.

**Data structures — the stream layer.** Shares the process-global `Stream` registry described
under `OpenRead`; the only difference from `OpenWrite` is the `fopen` mode, so
the returned object, the id handle, and the `Close`/`atexit` lifecycle are identical.

**Complexity / limits.** `O(1)`. Pure ANSI C99 (`fopen`), no POSIX guards. `ATTR_PROTECTED`.
