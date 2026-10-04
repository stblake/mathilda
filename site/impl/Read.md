---
source: src/io/read.c
---
**Algorithm.** `builtin_read` is the single-object reader that `ReadList` loops.
It peels trailing option rules (`RecordSeparators`, `WordSeparators`, `TokenWords`,
`NullRecords`, `NullWords`, defaulting from `Options[Read]`) with `options_extract`, then takes
the type spec — a default `Expression`, a single type, or a nested structure. The spec is
**validated before the file is touched** (`read_spec_valid` via `spec_scan`), so a malformed
spec never opens a stream. The first argument is resolved by `resolve_input`: an `InputStream`
object, or a `"file"`/`File["file"]` name that is auto-opened and **left open** (so successive
`Read["file", ...]` calls advance the same point; `Read::openx`/`noopen` → `$Failed`).

Reading runs over a `ReadCursor` into the stream's buffer. `read_structure` walks the spec
depth-first: a type leaf dispatches to `read_one`, and a compound spec (`{...}`, `Hold[...]`,
any head) rebuilds that head around its recursively-read children — so
`{{Number,Number},{Number,Number}}` reads a 2×2 matrix in one call. Leaf readers: `Byte`
(one byte as `0`–`255`), `Character` (one byte as a 1-char string), `Word`/`Record`/`String`
(separator-, record-, and line-delimited text, respectively), `Number`/`Real`
(`read_number`: the next word token is `parse_expression`'d and `evaluate`'d; a malformed token
prints `Read::readn` and yields `$Failed`; `Real` forces the value through `N`), and
`Expression` (`read_expression`: `parse_next_expression` reads one complete top-level form
**unevaluated**, and the assembled result is evaluated by the surrounding evaluator — so
`Read[s, Expression]` evaluates while `Read[s, Hold[Expression]]` stays held). The stream's
`pos` is advanced to the cursor afterwards.

**Data structures.** A `ReadCfg` (borrowed separator-string arrays from the options) and a
`ReadCursor` (`{buf, len, pos}`) over the input slot's resident buffer in the shared stream
registry (see `OpenRead`); the shared per-type readers (`next_word`,
`next_record`, `next_string_line`) and the longest-match separator scan (`sep_len_at`). A
`dup_n` copies substrings rather than using POSIX `strndup` (C99 portability, SPEC.md §10).

**Complexity / limits.** `O(bytes consumed)` per call. Past end of file a leaf returns
`EndOfFile`; within a partly-read structure the unread trailing slots are padded with
`EndOfFile`, and a read begun at EOF returns a bare `EndOfFile`. An invalid type spec leaves
the call unevaluated. Pure ANSI C99. `ATTR_PROTECTED`.
