---
source: src/io/read.c
---
**Definition.** `Expression` is a **type specification** used by `Read` and `ReadList`.
It has no builtin and no value — it is an inert token that tells the reader to parse one
complete Mathilda expression from the stream. It is the **default** type when none is
given (both `read.c` and `readlist.c` fall back to `SYM_Expression`).

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_Expression`). `read.c`'s
`type_from_symbol` maps it to the internal `RT_EXPRESSION` reader. The leaf is read
unevaluated and the assembled result is then evaluated, so `Read[s, Expression]`
evaluates the expression it read while `Read[s, Hold[Expression]]` keeps it in raw form.

**Usage & limits.** Meaningful only as a `Read`/`ReadList` type argument; evaluated on
its own it returns itself. With no explicit type, `ReadList["file"]` reads every
remaining expression. Past end of file the reader returns `EndOfFile`.
