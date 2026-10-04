---
source: src/io/read.c
---
**Definition.** `Number` is a type specification used in `Read` and `ReadList`: it reads a number, returned as an integer when the token has no decimal point or exponent and as an approximate number otherwise.
It is an inert symbol — no builtin, no value, no attributes — declared with its
docstring in `info.c` and recognised by the reader.

**Representation.** A bare `EXPR_SYMBOL`. The read-type scanner in `src/io/read.c`
maps it to an internal read-type tag (`if (n == SYM_Number) *out = RT_...`), after
which the matching reader pulls the value from the stream.  As a plain symbol
its `Head` is `Symbol` and it carries no attributes.

**Usage & limits.** Meaningful only as the type argument of `Read`/`ReadList` over a
file or an open stream; standing alone it is just a symbol. an integer-looking token reads back as an `Integer` and a decimal or exponent token as a `Real`.
