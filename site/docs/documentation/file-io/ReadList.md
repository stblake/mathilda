# ReadList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ReadList["file"]`**

reads all remaining expressions from a file and returns them as a list.

**`ReadList["file", type]`**

reads objects of the given type until end of file, returning the list read.

**`ReadList["file", {type1, type2, ...}]`**

reads one object of each type per pass, grouping each pass into a sublist, until end of file. If end of file is reached partway through a pass, the unread slots are filled with EndOfFile.

**`ReadList["file", types, n]`**

reads only the first n objects (or passes) of the specified types.

<details>
<summary>Notes</summary>

Types: Byte (integer code), Character (one-character string), Expression (a complete evaluated expression), Number (integer, or real if it has a decimal point), Real (always an approximate number, C/Fortran E notation accepted), Record (characters up to a record separator), String (a line), Word (characters delimited by word separators). Options RecordSeparators, WordSeparators, TokenWords, NullRecords and NullWords control tokenisation; see Options\[ReadList\]. ReadList also accepts an open InputStream. A named file that is not already open is opened and closed by ReadList. Returns $Failed if the file cannot be opened.

</details>

## Examples

_No verified examples yet for this function._

## Algorithm

readlist.c — ReadList[] file-reading builtin.

```text
ReadList["file"]                 read all remaining Expressions -> List
ReadList["file", type]           read `type` until EOF -> flat List
ReadList["file", {t1,...,tk}]    read one of each type per pass -> List of
                                 k-element sublists (EndOfFile-padded when
                                 EOF arrives partway through a pass)
ReadList["file", types, n]       stop after n top-level objects / passes
ReadList[InputStream[...], ...]  read from an already-open input stream
```

ReadList is a sequence of calls to the shared single-object reader: it loops

```text
read_structure() (read.c) to end of file, collecting each result.  The read
```

types, separator options, and EndOfFile padding all come from that one

```text
engine, which Read[] also uses.  A filename that is not already open is opened
```

and closed by ReadList (there is no user-visible stream in that case).

Options (trailing rules, defaults from Options[ReadList]): RecordSeparators, WordSeparators, TokenWords, NullRecords, NullWords.

ANSI C99 only; no POSIX symbols (SPEC.md §10).

## Implementation notes

- `Protected`.
- `$Failed` (with a `ReadList::noopen` diagnostic) if the file cannot be opened.
- A malformed `Number`/`Real` token prints `ReadList::readn`, contributes `$Failed` in that position, and reading continues.
- When end of file is reached partway through a `{type_1, ...}` pass, the unread slots of that final pass are filled with `EndOfFile`.
- `ReadList` is a sequence of `Read` calls: it loops the shared reading engine (`src/io/read.c`) to end of file. It also accepts an open `InputStream` (reading from its current point and leaving it open); a named file that is not already open is opened and closed by `ReadList`.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [Get](../../file-io/Get/), [Byte](../../other-advanced/Byte/), [Character](../../other-advanced/Character/), [Word](../../other-advanced/Word/), [Record](../../other-advanced/Record/), [String](../../other-advanced/String/), [Number](../../other-advanced/Number/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_readlist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_readlist.c)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)
