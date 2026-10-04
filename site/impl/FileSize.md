---
source: src/files.c
---
**Algorithm.** `builtin_filesize` accepts exactly one `EXPR_STRING` argument (anything else —
wrong arity, a symbol, a non-string atom — returns `NULL`, leaving the call unevaluated so
symbolic arguments flow through). It calls `stat()` on the path and, on success, returns
`st.st_size` as an `Integer` (`expr_new_integer((int64_t)st.st_size)`) — a plain byte count,
not a `Quantity`. `stat` (not `lstat`) follows symbolic links, so the size reported is the
target file's.

When the path cannot be `stat`'d (most often because nothing is there) it emits
`FileSize::nffil` through `fs_msg` and returns `$Failed`. `fs_msg` is the local diagnostic
funnel: it calls `mth_msg_note_fired()` (so `Check[]` sees the message) and returns early under
`mth_msg_suppressed()` (so `Quiet[]` silences it), then `vfprintf`s to `stderr` — the same
`Quiet`/`Check`-aware pattern as `dt_msg` in `src/datetime.c`.

**Data structures.** None beyond a `struct stat`; the path is read directly out of the argument
string. `_POSIX_C_SOURCE 200809L` is defined before any include so `stat` is visible under
glibc's `-std=c99` (SPEC.md §10).

**Complexity / limits.** `O(1)` — one `stat` syscall. The path is interpreted relative to the
current working directory; `$Path` is not searched. `ATTR_PROTECTED`.
