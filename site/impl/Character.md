---
source: src/io/read.c
---
**Definition.** `Character` is a **type specification** used by `Read` and `ReadList`.
It has no builtin and no value — it is an inert token that tells the reader to consume
exactly one character from the stream and return it as a one-character string.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_Character`). `read.c`'s
`type_from_symbol` maps it to the internal `RT_CHARACTER` reader, one of the eight read
types (`Byte`, `Character`, `Expression`, `Number`, `Real`, `Record`, `String`,
`Word`). It may appear on its own or nested inside a type-structure argument
(`{Character, Character}`, `Hold[...]`, any head), which the reader fills depth-first.

**Usage & limits.** Meaningful only as a `Read`/`ReadList` type argument; evaluated on
its own it simply returns itself. Past end of file the reader returns `EndOfFile`.
Contrast `Byte`, which returns the raw integer code 0-255 rather than a string.
