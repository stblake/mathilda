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

## Examples (16)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (11)

A three-line data file

```mathematica
In[1]:= Put[10, 20, 30, "/tmp/mathilda_io_read.txt"]
```

```mathematica
In[2]:= str = OpenRead["/tmp/mathilda_io_read.txt"]
Out[2]= InputStream["/tmp/mathilda_io_read.txt", 1]

In[3]:= Read[str, Number]
Out[3]= 10
```

Read the next two numbers into a list

```mathematica
In[4]:= Read[str, {Number, Number}]
Out[4]= {20, 30}
```

Past the end of the file

```mathematica
In[5]:= Read[str, Number]
Out[5]= EndOfFile
```

```mathematica
In[6]:= Close[str]
Out[6]= "/tmp/mathilda_io_read.txt"

In[7]:= str = OpenWrite["/tmp/mathilda_io_read2.txt"]; WriteString[str, "1 2\n3 4\nx + 1\n"]; Close[str];

In[8]:= ins = OpenRead["/tmp/mathilda_io_read2.txt"];
```

A 2x2 matrix in one call

```mathematica
In[9]:= Read[ins, {{Number, Number}, {Number, Number}}]
Out[9]= {{1, 2}, {3, 4}}
```

Read the next expression without evaluating it

```mathematica
In[10]:= Read[ins, Hold[Expression]]
Out[10]= Hold[x + 1]
```

```mathematica
In[11]:= Close[ins]
Out[11]= "/tmp/mathilda_io_read2.txt"
```

### Applications (5)

A small data file

```mathematica
In[12]:= str = OpenWrite["/tmp/mathilda_read.txt"]; WriteString[str, "1 2\n3 4\nx + 1\n"]; Close[str];
```

```mathematica
In[13]:= ins = OpenRead["/tmp/mathilda_read.txt"];
```

Read a 2x2 matrix in one call

```mathematica
In[14]:= Read[ins, {{Number, Number}, {Number, Number}}]
Out[14]= {{1, 2}, {3, 4}}
```

Read the next expression WITHOUT evaluating it

```mathematica
In[15]:= Read[ins, Hold[Expression]]
Out[15]= Hold[x + 1]
```

```mathematica
In[16]:= Close[ins]
Out[16]= "/tmp/mathilda_read.txt"
```

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

- `Protected`.
- Returns `EndOfFile` once past end of file. Within a partly-read `{type_1, ...}` structure, the unread trailing slots are filled with `EndOfFile`; a read that begins already at end of file returns a bare `EndOfFile`.
- Returns `$Failed` for a stream that is not open, or (with `Read::readn`) for a token that is not of the requested numeric type.
- The `Expression` leaf is read **unevaluated** and the constructed result is then evaluated by the surrounding evaluator — so `Read[s, Expression]` evaluates while `Read[s, Hold[Expression]]` stays held.

**Attributes:** `Protected`.

## References

**See also:** [ReadList](../../file-io/ReadList/), [InputStream](../../other-advanced/InputStream/), [OpenRead](../../file-io/OpenRead/), [RecordSeparators](../../other-advanced/RecordSeparators/), [WordSeparators](../../other-advanced/WordSeparators/), [TokenWords](../../other-advanced/TokenWords/), [NullRecords](../../other-advanced/NullRecords/), [NullWords](../../other-advanced/NullWords/)

- Source: [`src/io/read.c`](https://github.com/stblake/mathilda/blob/main/src/io/read.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)

## Notes & additional examples

### Notes

`Read` reads **one** object (or one nested type structure) from a stream and
advances its current point, so successive `Read` calls walk the file;
`ReadList` is this primitive looped to end of file. The type may be
`Byte`, `Character`, `Word`, `Record`, `String`, `Number`, `Real`, or
`Expression`, and a *structure* of types (`{{Number, Number}, ...}`, `Hold[...]`,
any head) is filled depth-first — the matrix above is read with a single call.

The `Expression` leaf is read unevaluated and the assembled result is then
evaluated, so `Read[s, Expression]` evaluates while `Read[s, Hold[Expression]]`
keeps the raw form. Past end of file `Read` returns `EndOfFile`; trailing slots of
a partly-read structure are padded with `EndOfFile`. The separator options
(`RecordSeparators`, `WordSeparators`, `TokenWords`, `NullRecords`, `NullWords`)
default from `Options[Read]`, so `SetOptions[Read, ...]` works.
