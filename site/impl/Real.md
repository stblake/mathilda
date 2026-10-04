---
source: src/io/read.c
---
**Definition.** `Real` has two inert roles, both without a builtin. First, it is the
**head of approximate real numbers**: `Head[3.14]` is `Real`, so `Real` names the type
of machine- and arbitrary-precision floating-point values (matched by `_Real` in
patterns). Second, it is a **type specification** for `Read` and `ReadList`, telling the
reader to read a number and always return it as an approximate (inexact) number, even
when the token has no decimal point.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_Real`). As a value head it is
attached to every `EXPR_REAL` / `EXPR_MPFR` node (the printer shows these as `3.14`,
not `Real[...]`). As a read type, `read.c`'s `type_from_symbol` maps it to the internal
`RT_REAL` reader, which accepts C/Fortran `E`-notation. The related `Number` type reads
an integer when the token is integral and a real otherwise; `Real` always coerces to a
real.

**Usage & limits.** As a read type it is meaningful only inside a `Read`/`ReadList`
call; on its own `Real` evaluates to itself. Past end of file the reader returns
`EndOfFile`.
