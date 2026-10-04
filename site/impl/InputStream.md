---
source: src/io/streams.c
---
**Definition.** `InputStream["name", n]` is the object that represents an open input
stream, as returned by `OpenRead`. It has no builtin and no value of its own — it is an
inert stream handle. `"name"` is the source (a file name), and the integer `n` is an
internal handle into the stream registry.

**Representation.** A two-argument `EXPR_FUNCTION` with head `InputStream` (interned
`SYM_InputStream`). `src/io/streams.c` builds it from an entry in the open-stream
registry (`InputStream[name, id]` for an input stream, `OutputStream[name, id]` for an
output stream) and resolves it back to that registry entry whenever `Read`, `ReadList`,
`SetStreamPosition`, `Skip` or `Close` is handed a stream. `Read`/`ReadList` also accept
a bare `"file"` or `File["file"]`, which they open, read, and leave open as an
`InputStream`.

**Usage & limits.** The handle is valid only while the underlying stream is open;
`Close` invalidates it. Reading from an already-open `InputStream` continues from its
current point and leaves it open (unlike a named-file read, which `ReadList` opens and
closes itself). The companion output object is `OutputStream`.
