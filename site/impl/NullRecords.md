---
source: src/io/read.c
---
**Definition.** `NullRecords` is an **option** for `Read` and `ReadList` specifying
whether a null (empty) record is assumed between two adjacent record separators. It has
no builtin and no value of its own — an inert, `Protected` keyword carried in the
reader's option list, with default `False`.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_NullRecords`). Both `read.c` and
`readlist.c` seed their default option set with `rule_of(SYM_NullRecords, False)`. When
`True`, the `Record` reader keeps the empty field that lies between two consecutive
`RecordSeparators`, so `"a,b,,c"` split on `","` yields `{"a", "b", "", "c"}` rather
than `{"a", "b", "c"}`.

**Usage & limits.** `Protected`; meaningful only as a `Read`/`ReadList` option in
combination with the `Record` read type and `RecordSeparators`. It has a registered
default, so `SetOptions` on it works. The word-level analogue is `NullWords`.
