---
source: src/io/streams.c
---
**Definition.** `OutputStream["name", n]` is the inert object representing an open
output stream, as returned by `OpenWrite` or `OpenAppend`. The string is the
stream's name (a file name, or `"stdout"`/`"stderr"`); the integer `n` is an
internal handle into the process-global stream registry. It is `Protected`, has no
builtin of its own, and its docstring lives in `info.c`.

**Representation.** `src/io/streams.c` keeps a process-global registry of open
streams. Each open stream is handed to the language as an `OutputStream[name, id]`
object (its input counterpart is `InputStream[name, id]`), where `id` addresses the
registry slot holding the live `FILE*`. The write and stream builtins (`Write`,
`WriteString`, `Close`, ...) resolve a target by looking that `id` up in the
registry, matching the head against the interned `SYM_OutputStream`. As a
two-argument expression its `Head` is `OutputStream`.

**Usage & limits.** It is a handle, not data: you pass it to the write and stream
builtins, and `Close` frees its registry slot. The integer handle is an internal
detail, meaningful only to the registry within the same session. Constructing one by
hand does **not** open a stream — a live object comes from `OpenWrite`/`OpenAppend`.
