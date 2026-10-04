---
source: src/io/read.c
---
**Definition.** `EndOfFile` is the **sentinel symbol** that `Read` returns when a stream
has no more to read. It is a bare symbol, not a function — no builtin, no DownValues, not
even the `Protected` attribute — carrying only the docstring set in `info_init`
(`src/info.c`). It is produced, not consumed, by the reader `src/io/read.c`.

**Representation.** `EndOfFile` stays an inert `EXPR_SYMBOL`; evaluated on its own it is
just `EndOfFile`, and `Head[EndOfFile]` is `Symbol`. The reader constructs it with
`expr_new_symbol(SYM_EndOfFile)` at the two points where a read runs out of input (the
scalar case, and the fallthrough that ends a `ReadList` loop).

**Usage & limits.** When `Read[stream, type]` is called at or past end of file it returns
`EndOfFile`; `ReadList` uses the same sentinel to fill the unread trailing slots of a type
sequence that end-of-file truncated, and stops its loop when it is produced. Code reading a
stream incrementally tests for it to know when to stop — `While[(x = Read[s, Word]) =!=
EndOfFile, ...]`. It is a plain marker value with no numeric or structural meaning of its
own.
