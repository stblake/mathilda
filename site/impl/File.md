---
source: src/io/streams.c
---
**Definition.** `File["name"]` is a **symbolic wrapper** for a file name, accepted wherever
`Read`, `ReadList`, `OpenRead`, `OpenWrite`, `OpenAppend` and `Close` take a file. It is an
inert head, not a function: `streams_init` in `src/io/streams.c` marks `File` `Protected`
(the docstring is set centrally in `src/info.c`), and there is no builtin that rewrites it —
`File["x"]` evaluates to itself.

**Representation.** `File["name"]` is an ordinary `EXPR_FUNCTION` with head `SYM_File` and
a single string argument, so `Head[File["data.txt"]]` is `File` and
`FullForm[File["data.txt"]]` is `File["data.txt"]`. The I/O layer never stores a `File`
object itself; the helper `stream_filename_arg` (`src/io/streams.c`) unwraps it, accepting
either a bare string or a one-argument `File[...]` whose argument is a string, and returns
the underlying path.

**Usage & limits.** `File` lets a path be passed as an explicit, typed object rather than a
raw string — `ReadList[File["data.txt"], Word]` behaves exactly like
`ReadList["data.txt", Word]`. The argument must be a single string; other forms are not
unwrapped. It is a name wrapper only: it does not open, test or resolve the file, and holds
no handle (an open stream is an `InputStream`/`OutputStream` object instead).
