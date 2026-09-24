# Read

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Read[stream]`**

reads one expression from an open input stream and returns it.

**`Read[stream, type]`**

reads one object of the specified type.

**`Read[stream, {type1, type2, ...}]`**

reads a sequence of objects of the specified types into a list.

**`Read[stream, structure]`**

reads into any nested type structure (for example {{Number, Number}, ...} or Hold\[Expression\]), filled by a depth-first traversal.

<details>
<summary>Notes</summary>

The stream may be an InputStream from OpenRead, a "file", or File\["file"\]; a named file that is not already open is opened and left open, so successive Read calls advance the same current point. Types are Byte, Character, Expression, Number, Real, Record, String and Word, as for ReadList; options RecordSeparators, WordSeparators, TokenWords, NullRecords and NullWords control tokenisation. Read returns EndOfFile past end of file, and $Failed for a stream that is not open or a token that is not of the requested type.

</details>

## Examples

_No verified examples yet for this function._

## Algorithm

read.c — the shared file-reading engine and the Read[] builtin.

```text
Read[stream]                 read one Expression from a stream
Read[stream, type]           read one object of the given type
Read[stream, {t1, t2, ...}]  read one object of each type into a list
Read[stream, structure]      any nested type-structure (Hold[...], f[...],
                             {{Number,Number},...}); filled depth-first

`stream` is an InputStream[...] object (from OpenRead), a bare "file" string,
or File["file"].  A filename that is not already open is opened and left open,
```

so successive Read calls advance the same current point (Close finishes it).

Read returns EndOfFile once past end of file, and $Failed for a malformed

```text
numeric token or a stream that is not open.  The per-type readers and the
```

separator handling are shared with ReadList (readlist.c), which loops read_structure() to end of file.

ANSI C99 only (fopen/fread via streams.c, strtod via the parser); no POSIX symbols, so no feature-test guards are needed (SPEC.md §10).

## Implementation notes

- `Protected`.
- Returns `EndOfFile` once past end of file. Within a partly-read `{type_1, ...}` structure, the unread trailing slots are filled with `EndOfFile`; a read that begins already at end of file returns a bare `EndOfFile`.
- Returns `$Failed` for a stream that is not open, or (with `Read::readn`) for a token that is not of the requested numeric type.
- The `Expression` leaf is read **unevaluated** and the constructed result is then evaluated by the surrounding evaluator — so `Read[s, Expression]` evaluates while `Read[s, Hold[Expression]]` stays held.

**Attributes:** `Protected`.

## References

**See also:** [ReadList](../../file-io/ReadList/), [InputStream](../../other-advanced/InputStream/), [OpenRead](../../file-io/OpenRead/), [RecordSeparators](../../other-advanced/RecordSeparators/), [WordSeparators](../../other-advanced/WordSeparators/), [TokenWords](../../other-advanced/TokenWords/), [NullRecords](../../other-advanced/NullRecords/), [NullWords](../../other-advanced/NullWords/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)
