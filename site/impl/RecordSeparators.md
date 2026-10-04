---
source: src/io/read.c
---
**Definition.** `RecordSeparators` is an **option** for `Read` and `ReadList` giving the
list of strings that delimit records (for the `Record` read type). It has no builtin and
no value of its own — it is an inert, `Protected` keyword carried in the reader's option
list, with default `{"\r\n", "\n", "\r"}`.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_RecordSeparators`). Both
`read.c` and `readlist.c` seed their default option set with
`rule_of(SYM_RecordSeparators, {"\r\n", "\n", "\r"})`; a user `RecordSeparators ->
{...}` rule overrides it. The reader's `next_record` scanner splits on any of these
strings. The companion `NullRecords` option decides whether an empty field between two
adjacent separators is kept.

**Usage & limits.** `Protected`; meaningful only as a `Read`/`ReadList` option. Because
it has a registered default, `SetOptions[Read, RecordSeparators -> ...]` works. Related
tokenisation options are `WordSeparators`, `TokenWords`, `NullRecords` and `NullWords`.
