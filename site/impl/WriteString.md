---
source: src/io/streams.c
---
**Algorithm.** `builtin_writestring` resolves its first argument with `resolve_output` exactly
as `Write` does (an `OutputStream`, or a name auto-opened truncating and left open;
`$Failed` otherwise). It then writes each remaining argument **verbatim**: a string is `fputs`'d
raw, with no surrounding quotes, and a non-string argument is rendered with `expr_to_string`
first. It adds **no separators between arguments and no trailing newline**, then `fflush`es and
returns `Null`.

**Data structures — the stream layer.** Writes to the `FILE*` of an output slot in the
process-global `Stream` registry described under `OpenRead`. It is the exact-byte
counterpart to `Write`: `Write` prints input form and appends `'\n'`, whereas
`WriteString` emits precisely the bytes given — so a caller lays out CSV lines, headers, or any
newline structure explicitly (e.g. `WriteString[s, "1 2 3\n4 5 6\n"]`).

**Complexity / limits.** `O(total output length)`; `fflush` per call. Pure ANSI C99
(`fputs`/`fflush`). `ATTR_PROTECTED`.
