---
source: src/core.c
---
**Definition.** `String` has two roles. Primarily it is the **head of string
objects**: a string leaf is an `EXPR_STRING`, and `Head` of one returns the symbol
`String`, so `_String` is the pattern that matches any string. Secondarily it is a
type specification in `Read` and `ReadList`, where it reads one line (up to a
newline). The symbol itself is inert — no builtin, no value; its docstring is in
`info.c`.

**Representation.** Strings are stored as `EXPR_STRING` nodes carrying a C string,
not as `String[...]` compounds; `String` is the *head* those leaves report, supplied
by the `Head` builtin in `src/core.c` (the `EXPR_STRING` case), which is why
`StringQ`/`MatchQ[..., _String]` work. As a read type, the scanner in `src/io/read.c`
maps `n == SYM_String` to the `RT_STRING` read tag, whose reader returns the next
line as a string. The bare symbol `String` is itself an `EXPR_SYMBOL`, so
`Head[String]` is `Symbol`.

**Usage & limits.** The head role is pervasive — any pattern, type test, or
`Cases`/`Select` over strings leans on it. The read-type role needs a file or open
stream, so it cannot be exercised on a bare literal. `String` performs no
computation on its own in either role.
