---
source: src/io/read.c
---
**Definition.** `Byte` is a type specification used in `Read` and `ReadList`: it reads a single byte, returned as an integer code from 0 to 255.
It is an inert symbol — no builtin, no value, no attributes — declared with its
docstring in `info.c` and recognised by the reader.

**Representation.** A bare `EXPR_SYMBOL`. The read-type scanner in `src/io/read.c`
maps it to an internal read-type tag (`if (n == SYM_Byte) *out = RT_...`), after
which the matching reader pulls the value from the stream. Because `Byte` reads the raw stream, it returns the numeric code of each octet, not a parsed number. As a plain symbol
its `Head` is `Symbol` and it carries no attributes.

**Usage & limits.** Meaningful only as the type argument of `Read`/`ReadList` over a
file or an open stream; standing alone it is just a symbol. it returns raw byte codes (so the character "A" reads as 65), unlike `Number`, which parses numeric tokens.
