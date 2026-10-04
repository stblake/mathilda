---
source: src/io/read.c
---
**Definition.** `Record` is a type specification used in `Read` and `ReadList`: it reads a sequence of characters delimited by record separators (see `RecordSeparators`), returned as a string.
It is an inert symbol — no builtin, no value, no attributes — declared with its
docstring in `info.c` and recognised by the reader.

**Representation.** A bare `EXPR_SYMBOL`. The read-type scanner in `src/io/read.c`
maps it to an internal read-type tag (`if (n == SYM_Record) *out = RT_...`), after
which the matching reader pulls the value from the stream. The record boundaries come from the `RecordSeparators` option, not from `Byte`/`Number` parsing. As a plain symbol
its `Head` is `Symbol` and it carries no attributes.

**Usage & limits.** Meaningful only as the type argument of `Read`/`ReadList` over a
file or an open stream; standing alone it is just a symbol. each record is returned as a `String`, split on the `RecordSeparators` strings.
