---
source: src/io/readlist.c
---
**Algorithm.** `builtin_readlist` is `Read`'s single-object engine **looped to end
of file**. It peels the same trailing option rules (`RecordSeparators`, `WordSeparators`,
`TokenWords`, `NullRecords`, `NullWords`, from `Options[ReadList]`) with `options_extract`,
validates the type spec with `read_spec_valid` before opening anything, and reads an optional
non-negative count `n`. The source is resolved by `resolve_input`: an already-open
`InputStream` is read from its current point and **left open**, while a `"file"`/`File["file"]`
name is opened *and closed by ReadList itself* (`opened_here`), so no user-visible stream is
created in that case (`ReadList::openx`/`noopen` → `$Failed`).

The loop repeatedly calls the shared `read_structure` over a `ReadCursor`, appending each
result to a `malloc`'d array (doubled as needed) until end of file (a `NULL` at a pass
boundary), until the count `n` is reached, or — defensively — until a pass makes no progress.
The collected array is wrapped in `List[...]`. With a type list `{t1, ..., tk}` each pass reads
one object of each type into a `k`-element sublist, and `EndOfFile` pads a pass that hits EOF
partway through.

**Data structures.** A `ReadCfg` + `ReadCursor` over the stream's resident buffer (shared with
`Read`; see `OpenRead` for the registry), and a grown `Expr**` collection array.
All read types, separator handling, and `EndOfFile` padding come from the one engine in
`read.c`.

**Complexity / limits.** `O(file size)` time; the result list holds one entry per object/pass,
so memory is `O(objects read)`. The distinction from `Get` is that results are
**collected** (not just evaluated for effect), and the non-expression read types let a data
file be parsed field by field. Pure ANSI C99. `ATTR_PROTECTED`.
