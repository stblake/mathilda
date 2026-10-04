---
source: src/io/read.c
---
**Definition.** `NullWords` is an **option name** for `Read` and `ReadList` specifying
whether a null (empty) word is assumed between two adjacent word separators. It is a bare
option symbol, not a function — no builtin, no DownValues, not even the `Protected`
attribute — carrying only the docstring in `info_init` (`src/info.c`). The reader
`src/io/read.c` reads it out of the option list (its default `NullWords -> False` is part
of `Options[Read]` / `Options[ReadList]`).

**Representation.** `NullWords` stays an inert `EXPR_SYMBOL`, appearing only on the left of
a rule; `FullForm[NullWords -> True]` is `Rule[NullWords, True]`. The value is a boolean;
the default is `False`.

**Usage & limits.** With `NullWords -> False` (the default) the word reader skips an empty
field between adjacent separators — two spaces read as one gap. With `NullWords -> True` it
emits an empty string `""` for that field, so `"a  b"` reads as `{"a", "", "b"}`. It is one
of the tokenisation options alongside `WordSeparators`, `TokenWords`, `RecordSeparators`
and `NullRecords`, all defaulting from `Options[Read]` so `SetOptions[Read, ...]` adjusts
them globally.
