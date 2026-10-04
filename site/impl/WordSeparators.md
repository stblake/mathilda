---
source: src/io/read.c
---
**Definition.** `WordSeparators` is an option for `Read` and `ReadList` giving the
list of strings that separate words when the type specification is `Word`. The
default is `{" ", "\t"}`. It is an inert option keyword — no builtin and no value of
its own (it carries no attributes) — declared with its docstring in `info.c` and it
appears in `Options[Read]` / `Options[ReadList]` with its default.

**Representation.** A bare `EXPR_SYMBOL`. During a read, the options are gathered
into a small table (`src/io/read.c`, mirrored in `readlist.c`): the `"WordSeparators"`
row binds a local pointer `o_word` from the user's rules, which `read_cfg_build`
folds into the read configuration. The tokeniser then uses that separator list when
cutting `Word` tokens out of the stream (`RT_WORD`), falling back to the string-list
default when the option is not given.

**Usage & limits.** Meaningful only in a `Read`/`ReadList` call that reads `Word`
tokens; it does not affect `Number`, `Record` (which uses `RecordSeparators`) or
`String` (whole line) reads. The symbol performs no computation — it names the
option whose value the reader consults.
