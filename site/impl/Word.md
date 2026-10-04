---
source: src/io/read.c
---
**Definition.** `Word` is a **read type** token for `Read` and `ReadList`: it asks the
reader to return the next run of characters delimited by word separators as a string. It is
a bare symbol, not a function — no builtin, no DownValues, not even the `Protected`
attribute — carrying only the docstring set in `info_init` (`src/info.c`). Its meaning
lives in the reader `src/io/read.c`, where `read_type_of` maps the symbol `SYM_Word` to the
internal `RT_WORD` case.

**Representation.** `Word` stays an inert `EXPR_SYMBOL`; evaluated on its own it is just
`Word`, and `Head[Word]` is `Symbol`. It appears only as a type argument, e.g.
`Read[stream, Word]` or `ReadList["file", Word]`, and inside the type *structure* of a
nested read such as `ReadList["file", {Word, Number}]`.

**Usage & limits.** When the reader hits a `Word` slot it skips leading word separators
(by default space and tab, configurable through `WordSeparators`), then collects characters
up to the next separator and returns them as a `String`. `TokenWords` can force specific
strings to read as standalone words, and `NullWords -> True` keeps empty fields between
adjacent separators. A `Word` read past end of file returns `EndOfFile`. `Word` is distinct
from `Record` (a whole record up to a record separator) and `String` (a line).
