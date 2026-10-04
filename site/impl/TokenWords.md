---
source: src/io/read.c
---
**Definition.** `TokenWords` is an **option** for `Read` and `ReadList` giving a list of
strings that are to be read as separate words even when they are not surrounded by word
separators. It has no builtin and no value of its own — an inert, `Protected` keyword
carried in the reader's option list, with default `{}`.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_TokenWords`). Both `read.c` and
`readlist.c` seed their default option set with `rule_of(SYM_TokenWords, {})`. The
reader's `next_word` scanner checks, at each position, whether a `TokenWords` string
begins there (`token_len`); if so it emits that token as its own field and advances past
it, so an operator glued to its operands (`a+b` with `TokenWords -> {"+"}`) splits into
`a`, `+`, `b`.

**Usage & limits.** `Protected`; meaningful only as a `Read`/`ReadList` option for the
`Word` read type. It has a registered default, so `SetOptions[Read, TokenWords -> ...]`
works. Related tokenisation options are `WordSeparators`, `RecordSeparators`,
`NullWords` and `NullRecords`.
