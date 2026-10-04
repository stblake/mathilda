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

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (7)

Put writes one expression per line

```mathematica
In[1]:= Put[2, 3, 5, 7, 11, "/tmp/mathilda_io_rl.txt"]
```

Every number into a flat list

```mathematica
In[2]:= ReadList["/tmp/mathilda_io_rl.txt", Number]
Out[2]= {2, 3, 5, 7, 11}
```

No type: every remaining expression

```mathematica
In[3]:= ReadList["/tmp/mathilda_io_rl.txt"]
Out[3]= {2, 3, 5, 7, 11}
```

```mathematica
In[4]:= str = OpenWrite["/tmp/mathilda_io_rl2.txt"]; WriteString[str, "a 1\nb 2\nc 3\n"]; Close[str];
```

One of each type per pass

```mathematica
In[5]:= ReadList["/tmp/mathilda_io_rl2.txt", {Word, Number}]
Out[5]= {{"a", 1}, {"b", 2}, {"c", 3}}
```

```mathematica
In[6]:= str = OpenWrite["/tmp/mathilda_io_csv.txt"]; WriteString[str, "a,b,c\n"]; Close[str];
```

Comma-separated fields

```mathematica
In[7]:= ReadList["/tmp/mathilda_io_csv.txt", Word, WordSeparators -> {","}]
Out[7]= {"a", "b", "c"}
```

### Applications (6)

Put writes one expression per line

```mathematica
In[8]:= Put[2, 3, 5, 7, 11, "/tmp/mathilda_rl.txt"];
```

Read every number into a flat list

```mathematica
In[9]:= ReadList["/tmp/mathilda_rl.txt", Number]
Out[9]= {2, 3, 5, 7, 11}
```

```mathematica
In[10]:= str = OpenWrite["/tmp/mathilda_rl2.txt"]; WriteString[str, "a 1\nb 2\nc 3\n"]; Close[str];
```

One of each type per pass -> sublists

```mathematica
In[11]:= ReadList["/tmp/mathilda_rl2.txt", {Word, Number}]
Out[11]= {{"a", 1}, {"b", 2}, {"c", 3}}
```

```mathematica
In[12]:= str = OpenWrite["/tmp/mathilda_rl3.txt"]; WriteString[str, "a,b,c\n1,2,3\n"]; Close[str];
```

Comma-separated fields

```mathematica
In[13]:= ReadList["/tmp/mathilda_rl3.txt", Word, WordSeparators -> {","}]
Out[13]= {"a", "b", "c", "1", "2", "3"}
```

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

- `Protected`.
- `$Failed` (with a `ReadList::noopen` diagnostic) if the file cannot be opened.
- A malformed `Number`/`Real` token prints `ReadList::readn`, contributes `$Failed` in that position, and reading continues.
- When end of file is reached partway through a `{type_1, ...}` pass, the unread slots of that final pass are filled with `EndOfFile`.
- `ReadList` is a sequence of `Read` calls: it loops the shared reading engine (`src/io/read.c`) to end of file. It also accepts an open `InputStream` (reading from its current point and leaving it open); a named file that is not already open is opened and closed by `ReadList`.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [Get](../../file-io/Get/), [Byte](../../other-advanced/Byte/), [Character](../../other-advanced/Character/), [Word](../../other-advanced/Word/), [Record](../../other-advanced/Record/), [String](../../other-advanced/String/), [Number](../../other-advanced/Number/)

- Source: [`src/io/readlist.c`](https://github.com/stblake/mathilda/blob/main/src/io/readlist.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_readlist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_readlist.c)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)

## Notes & additional examples

### Notes

`ReadList["file", type]` reads objects of `type` until end of file into a flat
list; `ReadList["file", {t1, ..., tk}]` reads one of each type per pass and groups
each pass into a `k`-element sublist (handy for columnar data). With no type it
reads every remaining expression. A trailing integer `n` stops after `n` objects or
passes.

Unlike `Get`, which only evaluates each expression for its effect,
`ReadList` **collects** the results, and the non-expression read types let a data
file be parsed field by field. `ReadList` loops the shared reading engine that
`Read` uses, so the read types and the separator options
(`RecordSeparators`, `WordSeparators`, `TokenWords`, `NullRecords`, `NullWords`,
from `Options[ReadList]`) are exactly the same. A named file is opened and closed by
`ReadList`; an already-open `InputStream` is read from its current point and left
open.
