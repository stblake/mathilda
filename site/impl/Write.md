---
source: src/io/streams.c
---
**Algorithm.** `builtin_write` resolves its first argument with `resolve_output`: an
`OutputStream` object, or a `"file"`/`File["file"]` name that is **auto-opened** (truncating)
and left open; an `InputStream` or an unopenable name fails with `$Failed` (`Write::openx` /
`Write::noopen`). For each remaining argument it renders the expression with `expr_to_string`
(re-readable input form) and `fputs` it, then writes a single `'\n'` and `fflush`es. Returns
`Null`. Because its arguments are ordinary (evaluated) arguments, `1 + 1` is written as `2`;
`Hold[...]` is the way to write an unevaluated form.

**Data structures — the stream layer.** Writes to the `FILE*` of an output slot in the
process-global `Stream` registry described under `OpenRead`. The `fflush` after
each call is what makes the bytes immediately visible to `Read`/`ReadList`,
so a write-then-read round trip works without an intervening `Close`.

**Complexity / limits.** `O(total output length)`. The trailing newline per call is what makes
the output parse back one object per `Read`. Differs from `WriteString`,
which writes strings raw with no quoting and no added newline. Pure ANSI C99
(`fputs`/`fputc`/`fflush`). `ATTR_PROTECTED`.
