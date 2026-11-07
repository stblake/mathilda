# File I/O

Builtins implemented in `src/readwrite.c` (`Get`/`Put`/`PutAppend`), the `src/io/` stream layer — `src/io/read.c` (`Read`, the shared reading engine) and `src/io/streams.c` (`StringToStream`/`OpenRead`/`OpenWrite`/`OpenAppend`/`Close`/`Streams`/`StreamPosition`/`SetStreamPosition`/`Write`/`WriteString`, plus the inert `InputStream`/`OutputStream`/`File` objects) — `src/io/readlist.c` (`ReadList`), and `src/files.c` (`FileExistsQ`, `FileSize`, `FileExtension`, `FileBaseName`, `FileNameJoin`, `FileNameSplit`, `FilePrint`, `DirectoryName`, `FileNames`, `Directory`, `SetDirectory`, `ResetDirectory`, `DirectoryStack`, `Environment`, `ParentDirectory`, `FileNameTake`, `AbsoluteFileName`, `FileType`, `DirectoryQ`, `ToFileName`, `DeleteFile`, and the `$InputFileName` / `$HomeDirectory` / `$TemporaryDirectory` / `$BaseDirectory` / `$UserBaseDirectory` / `$PreferencesDirectory` / `$UserDocumentsDirectory` / `$PathnameSeparator` system variables).

## Get
Reads a sequence of Mathilda expressions from a file, evaluates each in order, and returns the value of the last one.
- `Get["filename"]`
- `Get["ctx`"]`: reads the file holding package `ctx`` (used by `Needs`).
- `<<file` / `<< "file"`: the shorthand form. The name is a file token, not an
  expression: letters, digits, `` !$*-./:?\_`~ `` and balanced `[...]`
  (`<<f[x]` is `Get["f[x]"]`), ended by whitespace, `;`, `,`, `)`, `]` or `}`;
  any other character right after it is a syntax error (`<<f+1`). `Get` binds
  tighter than any operator: `a<<f` is `Times[a, Get["f"]]`.

**Features**:
- `Protected`.
- **`$Path` search.** A relative file name that does not open as given is
  looked for under each directory of `$Path` in turn. A context name ``"A`B`"``
  is tried as `A/B.m`, `A/B.wl`, `A/B/init.m`, `A/B/init.wl`,
  `A/B/Kernel/init.m`, `A/B/Kernel/init.wl` under each `$Path` directory.
  Absolute paths are used as given.
- `$Path` is an ordinary assignable list of directory strings, default `{"."}`
  (`AppendTo[$Path, "dir"]`). Non-string and empty entries are skipped.
- Returns `$Failed`, with `Get::noopen` (Quiet/Check-aware), if the file cannot
  be found or opened.
- **Syntax errors.** At a syntax error `Get` emits `Syntax::sntx` (with the
  offending text and its line number), stops reading, and returns `Null`.
  The statements before the error keep their effect. The message goes through
  the funnel, so `Quiet` and `Check` see it, and the parser's own diagnostics
  are not printed as well. The statement-splitting rules, including
  `lhs := rhs;/; cond` on one line, are in `docs/spec/operators.md`.
- Used by the REPL bootstrap to load `src/internal/init.m` (and the rules it pulls in).
- Files conventionally end with `.m`.

## ReadList
Reads the objects contained in a file and returns them as a `List`, directed by a *read type* and modified by separator options. Unlike `Get`, the objects are **collected** (not just evaluated for their side effect), and non-expression types let a data file be read field by field.
- `ReadList["file"]` — read every remaining expression; equivalent to `ReadList["file", Expression]`.
- `ReadList["file", type]` — read objects of `type` until end of file → a flat list.
- `ReadList["file", {type_1, type_2, ...}]` — read one object of each type per pass, grouping each pass into a sublist, until end of file → a list of sublists.
- `ReadList["file", types, n]` — read only the first `n` objects (single type) or `n` passes (type list).

**Read types**:

| Type | Reads | Returns |
|------|-------|---------|
| `Byte` | one byte | integer code `0`–`255` |
| `Character` | one byte | a one-character string |
| `Word` | a run of non-separator characters | string |
| `Record` | characters up to a record separator | string |
| `String` | a line (up to a newline; a trailing `\r` is dropped) | string |
| `Number` | one word token | integer if it has no decimal point or exponent, otherwise a real |
| `Real` | one word token | always an approximate number |
| `Expression` | one complete top-level expression | the expression (evaluated, unless an enclosing structure such as `Hold` protects it) |

- `Number`/`Real` accept the C/Fortran scientific forms: `2.e5`, `2E5`, `1.5e-3`, and `2.*10^5` all denote `2×10^5`. `Real` always returns an approximate number; `Number` returns an integer only when the token has no explicit decimal point or exponent.

**Options** (trailing rules; the defaults are the single source of truth in `Options[ReadList]`, so `SetOptions[ReadList, ...]` works):
- `RecordSeparators -> {"\r\n", "\n", "\r"}` — the strings that terminate a `Record` (record separators also terminate a `Word`).
- `WordSeparators -> {" ", "\t"}` — the strings that separate `Word`/`Number`/`Real` tokens.
- `TokenWords -> {}` — strings always read as separate words even without surrounding separators (`TokenWords -> {"+"}` reads `a+b` as `{"a", "+", "b"}`).
- `NullRecords -> False` — when `True`, empty records between adjacent record separators are kept rather than skipped.
- `NullWords -> False` — when `True`, empty words between adjacent word separators are kept.

**Features**:
- `Protected`.
- `$Failed` (with a `ReadList::noopen` diagnostic) if the file cannot be opened.
- A malformed `Number`/`Real` token prints `ReadList::readn`, contributes `$Failed` in that position, and reading continues.
- When end of file is reached partway through a `{type_1, ...}` pass, the unread slots of that final pass are filled with `EndOfFile`.
- `ReadList` is a sequence of `Read` calls: it loops the shared reading engine (`src/io/read.c`) to end of file. It also accepts an open `InputStream` (reading from its current point and leaving it open); a named file that is not already open is opened and closed by `ReadList`.

#### Basic Examples

```mathematica
In[1]:= Put[2, 3, 5, 7, 11, "/tmp/mathilda_io_rl.txt"]  (* Put writes one expression per line *)

In[2]:= ReadList["/tmp/mathilda_io_rl.txt", Number]  (* every number into a flat list *)
Out[2]= {2, 3, 5, 7, 11}

In[3]:= ReadList["/tmp/mathilda_io_rl.txt"]  (* no type: every remaining expression *)
Out[3]= {2, 3, 5, 7, 11}

In[4]:= str = OpenWrite["/tmp/mathilda_io_rl2.txt"]; WriteString[str, "a 1\nb 2\nc 3\n"]; Close[str];

In[5]:= ReadList["/tmp/mathilda_io_rl2.txt", {Word, Number}]  (* one of each type per pass *)
Out[5]= {{"a", 1}, {"b", 2}, {"c", 3}}

In[6]:= str = OpenWrite["/tmp/mathilda_io_csv.txt"]; WriteString[str, "a,b,c\n"]; Close[str];

In[7]:= ReadList["/tmp/mathilda_io_csv.txt", Word, WordSeparators -> {","}]  (* comma-separated fields *)
Out[7]= {"a", "b", "c"}
```

## StringToStream
`StringToStream["string"]` opens an input stream that reads from a copy of the
string. The stream object is `InputStream["String", id]`.
- `Read`, `ReadList`, `ReadString`, `Skip`, `StreamPosition`,
  `SetStreamPosition` and `Close` work on it exactly as on a file stream.
  `Close` returns `"String"`.
- `Read[s, Hold[Expression]]` reads statements one at a time, unevaluated,
  which is how a loader can read a package header.
- A non-string argument emits `StringToStream::string` and gives `$Failed`.
- `Protected`.

## Skip
`Skip[stream, type]` reads and discards one object of the given type, exactly
as `Read` would read it, with the same separator options. `Skip[stream, type, n]`
skips `n` objects. It returns `Null`, or `EndOfFile` if the stream ends first.
`Protected`.

## ReadString
Reads text from a file or an input stream.
- `ReadString["file"]`: the complete contents of a file. The file is read whole
  each time, and no stream is registered.
- `ReadString[stream]`: reads from the stream's current point to the end, and
  advances the stream.
- `ReadString[src, "term"]` or `ReadString[src, {"t1", ...}]`: reads up to, but
  not including, the first occurrence of a terminator. The terminator is left
  unread.

**Features**:
- `Protected`.
- The bytes are returned exactly as they are; CRLF line ends stay `"\r\n"`.
- At end of file the result is `EndOfFile`.
- A file that cannot be opened emits `OpenRead::noopen` and gives `$Failed`.
- Implemented in `src/io/read.c`.

## Read
Reads **one** object (or one nested type structure) from an input stream, advancing the stream's current point, so successive `Read` calls return successive objects. `Read` is the primitive that `ReadList` loops.
- `Read[stream]` — read one expression; equivalent to `Read[stream, Expression]`.
- `Read[stream, type]` — read one object of the given type.
- `Read[stream, {type_1, type_2, ...}]` — read one object of each type into a list.
- `Read[stream, structure]` — read into any nested type structure, filled by a depth-first traversal. The structure may use any head, so `Read[stream, {{Number, Number}, {Number, Number}}]` reads a 2×2 matrix and `Read[stream, Hold[Expression]]` reads an expression without evaluating it.

`stream` is an `InputStream` object (from `OpenRead`), a `"file"` string, or `File["file"]`. A named file that is not already open is opened and **left open**, so successive `Read["file", ...]` calls advance the same current point (finish with `Close["file"]`).

The read types and the `RecordSeparators`/`WordSeparators`/`TokenWords`/`NullRecords`/`NullWords` options are exactly those of `ReadList` (defaults in `Options[Read]`).

**Features**:
- `Protected`.
- Returns `EndOfFile` once past end of file. Within a partly-read `{type_1, ...}` structure, the unread trailing slots are filled with `EndOfFile`; a read that begins already at end of file returns a bare `EndOfFile`.
- Returns `$Failed` for a stream that is not open, or (with `Read::readn`) for a token that is not of the requested numeric type.
- The `Expression` leaf is read **unevaluated** and the constructed result is then evaluated by the surrounding evaluator — so `Read[s, Expression]` evaluates while `Read[s, Hold[Expression]]` stays held.

#### Basic Examples

```mathematica
In[1]:= Put[10, 20, 30, "/tmp/mathilda_io_read.txt"]  (* a three-line data file *)

In[2]:= str = OpenRead["/tmp/mathilda_io_read.txt"]
Out[2]= InputStream["/tmp/mathilda_io_read.txt", 1]

In[3]:= Read[str, Number]
Out[3]= 10

In[4]:= Read[str, {Number, Number}]  (* read the next two numbers into a list *)
Out[4]= {20, 30}

In[5]:= Read[str, Number]  (* past the end of the file *)
Out[5]= EndOfFile

In[6]:= Close[str]
Out[6]= "/tmp/mathilda_io_read.txt"
```

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_io_read2.txt"]; WriteString[str, "1 2\n3 4\nx + 1\n"]; Close[str];

In[2]:= ins = OpenRead["/tmp/mathilda_io_read2.txt"];

In[3]:= Read[ins, {{Number, Number}, {Number, Number}}]  (* a 2x2 matrix in one call *)
Out[3]= {{1, 2}, {3, 4}}

In[4]:= Read[ins, Hold[Expression]]  (* read the next expression without evaluating it *)
Out[4]= Hold[x + 1]

In[5]:= Close[ins]
Out[5]= "/tmp/mathilda_io_read2.txt"
```

## OpenRead / OpenWrite / OpenAppend
Open a file and return a stream object addressing a persistent current point.
- `OpenRead["file"]` / `OpenRead[File["file"]]` — open for reading; returns `InputStream["file", n]`.
- `OpenWrite["file"]` — open for writing, truncating the file; returns `OutputStream["file", n]`.
- `OpenAppend["file"]` — open for writing at the end of the file; returns `OutputStream["file", n]`.

**Features**:
- `Protected`. Return `$Failed` (with an `Open*::noopen` diagnostic) if the file cannot be opened.
- The integer in the returned object is an internal handle into the stream registry; the object is inert and prints as itself.
- The registry is freed on `Close` or, for anything still open, at program exit (no leaks).
- Options: an output stream keeps `FormatType -> form` (the form `Write` renders in; `InputForm` by default, `FortranForm` and `CForm` supported) and `PageWidth -> n` (78 by default; a `FortranForm` or `CForm` stream breaks a wider expression into continuation lines, see `FortranForm`). Other options (`CharacterEncoding`, ...) are accepted and ignored.

#### Basic Examples

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_io_open.txt"]  (* OpenWrite truncates the file *)
Out[1]= OutputStream["/tmp/mathilda_io_open.txt", 1]

In[2]:= Write[str, x^2]

In[3]:= Close[str]
Out[3]= "/tmp/mathilda_io_open.txt"

In[4]:= app = OpenAppend["/tmp/mathilda_io_open.txt"]  (* OpenAppend keeps the prior contents *)
Out[4]= OutputStream["/tmp/mathilda_io_open.txt", 2]

In[5]:= Write[app, y^2]

In[6]:= Close[app]
Out[6]= "/tmp/mathilda_io_open.txt"

In[7]:= ins = OpenRead["/tmp/mathilda_io_open.txt"]
Out[7]= InputStream["/tmp/mathilda_io_open.txt", 3]

In[8]:= ReadList[ins]  (* both lines survive the append *)
Out[8]= {x^2, y^2}

In[9]:= Close[ins]
Out[9]= "/tmp/mathilda_io_open.txt"
```

## Close
Closes an open stream and returns its file name.
- `Close[stream]` — close an `InputStream`/`OutputStream` object.
- `Close["file"]` / `Close[File["file"]]` — close a stream opened for that file.

**Features**: `Protected`. Returns `$Failed` (with `Close::stream`) if the stream is not open.

#### Basic Examples

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_io_close.txt"]
Out[1]= OutputStream["/tmp/mathilda_io_close.txt", 1]

In[2]:= Streams["/tmp/mathilda_io_close.txt"]  (* the stream is open *)
Out[2]= {OutputStream["/tmp/mathilda_io_close.txt", 1]}

In[3]:= Close[str]  (* Close returns the file name *)
Out[3]= "/tmp/mathilda_io_close.txt"

In[4]:= Streams["/tmp/mathilda_io_close.txt"]  (* and removes it from the registry *)
Out[4]= {}
```

## Streams
Lists the currently open streams.
- `Streams[]` — all open `InputStream`/`OutputStream` objects.
- `Streams["file"]` — only the streams open for the named file.

**Features**: `Protected`.

#### Basic Examples

```mathematica
In[1]:= a = OpenWrite["/tmp/mathilda_io_a.txt"]
Out[1]= OutputStream["/tmp/mathilda_io_a.txt", 1]

In[2]:= b = OpenWrite["/tmp/mathilda_io_b.txt"]
Out[2]= OutputStream["/tmp/mathilda_io_b.txt", 2]

In[3]:= Streams[]  (* all open streams *)
Out[3]= {OutputStream["/tmp/mathilda_io_a.txt", 1], OutputStream["/tmp/mathilda_io_b.txt", 2]}

In[4]:= Streams["/tmp/mathilda_io_a.txt"]  (* just the streams for one file *)
Out[4]= {OutputStream["/tmp/mathilda_io_a.txt", 1]}

In[5]:= Close[a]
Out[5]= "/tmp/mathilda_io_a.txt"

In[6]:= Close[b]
Out[6]= "/tmp/mathilda_io_b.txt"
```

## StreamPosition / SetStreamPosition
Query and set the current point of a stream, as an integer byte offset.
- `StreamPosition[stream]` — the current position.
- `SetStreamPosition[stream, n]` — set the position to offset `n` (clamped into the file); returns the new position.
- `SetStreamPosition[stream, Infinity]` — move to the end of the stream.

**Features**: `Protected`. Return `$Failed` if the stream is not open.

#### Basic Examples

```mathematica
In[1]:= Put[100, 200, 300, "/tmp/mathilda_io_sp.txt"]  (* writes "100\n200\n300\n" *)

In[2]:= str = OpenRead["/tmp/mathilda_io_sp.txt"]
Out[2]= InputStream["/tmp/mathilda_io_sp.txt", 1]

In[3]:= Read[str, Number]
Out[3]= 100

In[4]:= StreamPosition[str]  (* byte offset after "100\n" *)
Out[4]= 4

In[5]:= SetStreamPosition[str, 0]  (* rewind to the start *)
Out[5]= 0

In[6]:= Read[str, Number]  (* the same first number again *)
Out[6]= 100

In[7]:= SetStreamPosition[str, Infinity]  (* jump to end of file: the byte length *)
Out[7]= 12

In[8]:= Close[str]
Out[8]= "/tmp/mathilda_io_sp.txt"
```

## Write / WriteString
Write to an output stream.
- `Write[stream, expr_1, expr_2, ...]` — write each expression in input form, then one newline.
- `WriteString[stream, str_1, ...]` — write the strings with no added quotes or newline; non-string arguments are written in input form.

`stream` may be an `OutputStream` object, a `"file"`, or `File["file"]`; a named file that is not already open is opened for writing (truncating) and left open.

The names `"stdout"` and `"stderr"` are the process's own output streams, not files: `WriteString["stdout", ...]` prints (alongside `Print`), `Write` there uses OutputForm (strings unquoted), and `Close` flushes them without closing.

**Features**:
- `Protected`. Return `$Failed` if the file cannot be opened.
- `Write` evaluates its expression arguments before writing (use `Hold[...]` to write an unevaluated form), and output is flushed after each call so it round-trips with `Read`/`ReadList`.

#### Basic Examples

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_io_w.txt"]
Out[1]= OutputStream["/tmp/mathilda_io_w.txt", 1]

In[2]:= Write[str, 1 + 1]  (* Write evaluates, prints input form, adds a newline *)

In[3]:= Write[str, a + b]

In[4]:= WriteString[str, "done\n"]  (* WriteString writes raw text, no added newline *)

In[5]:= Close[str]
Out[5]= "/tmp/mathilda_io_w.txt"

In[6]:= ReadList["/tmp/mathilda_io_w.txt"]  (* 1 + 1 was evaluated to 2 before writing *)
Out[6]= {2, a + b, done}
```

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_io_wh.txt"]; Write[str, Hold[2 + 2]]; Close[str];  (* Hold writes the unevaluated form *)

In[2]:= ReadList["/tmp/mathilda_io_wh.txt", String]
Out[2]= {"Hold[2 + 2]"}
```

## LoadModule
Loads an internal Mathilda source module, resolving its location independently of
the current working directory.
- `LoadModule["relpath"]` — `relpath` is relative to the source tree's
  `src/internal` directory (e.g. `"simp/FullSimplify.m"`).

**Features**:
- `Protected`. Returns `True` if the module was located and loaded (or had already
  been loaded), `False` otherwise.
- Resolution is independent of the current working directory and tries, in
  order: `$MATHILDA_HOME/<relpath>` (then `$MATHILDA_HOME/internal/<relpath>`, so
  pointing it at the directory *containing* `internal/` also works);
  `<exe_dir>/src/internal/<relpath>` and
  `<exe_dir>/../share/mathilda/internal/<relpath>` (relative to the running
  binary, so a relocated or installed executable still finds its modules);
  `$(PREFIX)/share/mathilda/internal/<relpath>` when built with a compile-time
  `MATHILDA_PREFIX`; and finally a CWD ladder (`src/internal/`,
  `../src/internal/`, `../../src/internal/`, `../../../src/internal/`). This
  works from the REPL (run at the repo root or from anywhere with the binary in
  place), from the test binaries (run from `tests/build/`), and from a binary
  copied to a `bin` directory with `MATHILDA_HOME` pointing at `src/internal`.
  The winning base directory is cached after the first successful lookup.
- Each module is loaded **at most once**, so repeated calls — and the lazy
  per-family loading used by [`FullSimplify`](simplification.md#fullsimplify) —
  never re-register rules.
- Generalises the bespoke fallback previously hard-coded for the CRC integral
  tables; `Get` (above) shares its file-reading core.

#### Basic Examples

```mathematica
In[1]:= LoadModule["simp/FullSimplify.m"]  (* load an internal module *)
Out[1]= True

In[2]:= LoadModule["simp/FullSimplify.m"]  (* already loaded: still True, not re-read *)
Out[2]= True

In[3]:= Quiet[LoadModule["no/such/module.m"]]  (* not found -> False *)
Out[3]= False
```

## Import
Reads a raster image file and returns an `Image`, or reads a text file.
- `Import["file"]` — decodes PNG, JPEG, BMP, GIF, TGA, PSD, HDR or PNM by content.
- `Import["file", "Image"]` — the same, stated explicitly.

**Text formats.** These also accept the list form `{"Text"}`. A `.txt` file
with no format given imports as `"Text"`.

| Format | Result for a file containing `"a\r\nb\n\n"` |
|---|---|
| `"String"` | the contents byte for byte: `"a\r\nb\n\n"` |
| `"Text"` | the contents with CRLF and lone CR read as `"\n"`, less one trailing line break: `"a\nb\n"` |
| `"Lines"` | a list of lines (split at LF, CRLF or a lone CR), with leading and trailing empty lines dropped: `{"a", "b"}` |

A file that cannot be read emits `Import::nffil` and gives `$Failed`.

**Features**:
- `Protected`.
- Samples are scaled by `1/255` into the unit interval, so the result is a `"Real"`
  image whatever the file's bit depth. The stored range of an image fixes what every
  downstream kernel's arithmetic means (see `ImageData`), and a range that depended on
  the file would make a Gaussian's scale depend on it too.
- The file's channel count is preserved: a grey file stays 1-channel and an RGBA file
  keeps its alpha. Forcing 3 channels would invent two for the first and silently
  discard transparency from the second.
- The result is packed and canonical — the same representation a filter produces, so an
  imported photograph needs no special-casing downstream.
- `$Failed` for a missing or malformed file. A path whose format is not handled at all
  stays unevaluated instead, which keeps `Import` from appearing to implement every
  format in existence.
- Decoding is by the vendored `stb_image` (public domain), so no system image library is
  a build requirement.

#### Basic Examples

```mathematica
In[1]:= img = Image[Table[{N[i/24], N[j/32], N[Mod[i + j, 8]]/8}, {i, 1, 24}, {j, 1, 32}], "Real"];

In[2]:= Export["/tmp/mathilda_doc.png", img]
Out[2]= /tmp/mathilda_doc.png

In[3]:= Import["/tmp/mathilda_doc.png"]
Out[3]= -Image-

In[4]:= ImageDimensions[Import["/tmp/mathilda_doc.png"]]
Out[4]= {32, 24}

In[5]:= ImageType[Import["/tmp/mathilda_doc.png"]]
Out[5]= Real
```

#### Scope

```mathematica
In[1]:= img = Image[Table[{N[i/24], N[j/32], N[Mod[i + j, 8]]/8}, {i, 1, 24}, {j, 1, 32}], "Real"];

In[2]:= ImageChannels[Import[Export["/tmp/mathilda_doc_g.png", Image[Table[N[i/16], {i, 1, 16}, {j, 1, 16}], "Real"]]]]
Out[1]= 1

In[2]:= ImageChannels[Import[Export["/tmp/mathilda_doc_a.png", Image[Table[{0.2, 0.4, 0.6, 0.8}, {i, 1, 8}, {j, 1, 8}], "Real"]]]]
Out[2]= 4

In[3]:= ImageDimensions[Import[Export["/tmp/mathilda_doc.jpg", img]]]
Out[3]= {32, 24}

In[4]:= ImageDimensions[Import[Export["/tmp/mathilda_doc.bmp", img]]]
Out[4]= {32, 24}

In[5]:= ImageDimensions[Import[Export["/tmp/mathilda_doc.tga", img]]]
Out[5]= {32, 24}

In[6]:= Import["/tmp/mathilda_doc_missing.png"]
Out[6]= $Failed
```

#### Properties & Relations

```mathematica
In[1]:= img = Image[Table[{N[i/24], N[j/32], N[Mod[i + j, 8]]/8}, {i, 1, 24}, {j, 1, 32}], "Real"];

In[2]:= (* PNG is lossless, so a round trip is exact to within half a quantisation level *)
Max[Abs[Flatten[ImageData[Import[Export["/tmp/mathilda_doc.png", img]]] - ImageData[img]]]] <= 1/510. + 1.*^-12
Out[1]= True

In[2]:= (* an imported image is packed, like every filter's result *)
Head[Part[Import["/tmp/mathilda_doc.png"], 1]]
Out[2]= NDArray

In[3]:= (* JPEG is lossy: the same round trip is bounded, not exact *)
0 < Mean[Flatten[Abs[ImageData[Import[Export["/tmp/mathilda_doc.jpg", img]]] - ImageData[img]]]] < 0.1
Out[3]= True

In[4]:= (* a format nothing here claims stays unevaluated, which is not the same failure as a missing file *)
Head[Import["/tmp/mathilda_doc.xyz"]]
Out[4]= Import
```

#### Applications

```mathematica
In[1]:= img = Image[Table[{N[i/24], N[j/32], N[Mod[i + j, 8]]/8}, {i, 1, 24}, {j, 1, 32}], "Real"];

In[2]:= (* filters compose with an imported image exactly as with a constructed one *)
ImageDimensions[GaussianFilter[Import["/tmp/mathilda_doc.png"], 2]]
Out[1]= {32, 24}

In[2]:= (* a pipeline written end to end: read, edge-detect, write *)
Export["/tmp/mathilda_doc_edges.png", EdgeDetect[Import["/tmp/mathilda_doc.png"]]]
Out[2]= /tmp/mathilda_doc_edges.png

In[3]:= ImageDimensions[Import["/tmp/mathilda_doc_edges.png"]]
Out[3]= {32, 24}
```

## Export
Writes an `Image` to a raster image file, or a `Graphics` object to an image file.
- `Export["file", image]` — the format comes from the file extension (PNG, JPEG, BMP, TGA).
- `Export["file", image, "PNG"]` — the format stated explicitly, for a name that does not
  carry one.
- `Export["file", plot]` where `plot` is the result of `Plot`, `ListPlot`, `Graphics`, …
  writes the graphic to **PDF, PNG, or JPEG** (again chosen from the extension, or stated
  as `"PDF"`/`"PNG"`/`"JPEG"`).

**Features**:
- `Protected`.
- Returns the file name, so `Import[Export[f, img]]` is a round trip that can be written
  as a single expression.

**Graphics export**:
- **PDF** is a resolution-independent **vector** file written by a small built-in PDF
  emitter — no external library, no display, so it works headless and in every build. It
  walks the graphics primitives directly (`Line`, `Point`, `Polygon`, `Disk`/`Circle`,
  `Rectangle`, `Arrow`, `Text`) with the `RGBColor`/`GrayLevel`/`Hue`/`CMYKColor`,
  `Opacity`, `Thickness`/`AbsoluteThickness`, `PointSize`/`AbsolutePointSize`,
  `Dashing`/`AbsoluteDashing` (`Dashed`, `Dotted`, `DotDashed`) and `Directive`
  directives, each `{...}` list scoping the directives set inside it, and draws a framed
  set of axes with "nice" ticks and numeric labels. It honours `PlotStyle` (the style a
  single-curve `Plot`/`ListPlot`/`ParametricPlot` is drawn in), `Prolog`/`Epilog`,
  `PlotLabel`, `AxesLabel` and `PlotLegends` (a swatch-and-label legend to the right of the
  plot). Coordinates may be exact (`1/2`, `Pi/4`, `Sqrt[2]`): they are converted the same
  way the on-screen renderer converts them. Text uses the PDF base-14 Helvetica, so no font
  is embedded. This is the recommended format for print and for the book.
- In the PDF, `Text[s, pos, {ox, oy}]` aligns as the reference does (`{-1, 0}` puts the left
  end of `s` at `pos`, `{0, 0}` centres it, using the Helvetica advance widths), and
  `Text[Style[s, n | FontSize -> n | colour, ...], ...]` sets that string's size and
  colour. `Arrowheads[s]` fixes the arrowhead length at `s` times the plot width, and an
  arrow's shaft stops inside its head rather than poking past the point.
  `AspectRatio -> Automatic` maps x and y with one scale (the page height follows the
  data unless `ImageSize -> {w, h}` fixes both, in which case the picture is centred).
- **PNG** and **JPEG** render through the graphics backend into an offscreen buffer, so the
  file is pixel-identical to the on-screen plot (the same axes, ticks, labels and text).
  They therefore need graphics support compiled in (`USE_GRAPHICS`) **and** a usable GUI
  session; with none (a headless box, `ssh`, cron) they return `$Failed` gracefully rather
  than crashing, while PDF still works. Resolution follows the `ImageSize` option: a width
  (default 800) with the height derived from the plot's `AspectRatio`, sized exactly as the
  on-screen window is, so an aspect-driven plot (`ArrayPlot`, `DensityPlot`, `ContourPlot`,
  …) fills its frame edge-to-edge instead of letterboxing inside a fixed canvas;
  `ImageSize -> {w, h}` pins both dimensions (then `AspectRatio` shapes the data inside that
  box). The pixels are encoded by the vendored `stb_image_write`, so JPEG output does not
  depend on which formats the Raylib build happens to support. Directive sizes mean the
  same thing in both formats: `Thickness`/`PointSize`/`Dashing` are fractions of the plot
  width and the `Absolute*` forms are printer's points (one pixel in a raster).
- A `Graphics3D` object (`Plot3D`, `ParametricPlot3D`, `ComplexPlot3D`, ...) exports to
  **PNG or JPEG** through the 3D renderer, with the same graphics-support/display
  requirement; it has no vector-PDF form (PDF of a 3D scene returns `$Failed`).
- Samples outside the unit interval are **clamped**, not wrapped. An unsharp mask
  legitimately overshoots and 8-bit output has nowhere to put the overshoot; wrapping
  would turn a bright highlight black, which reads as a bug in the filter rather than in
  the writer. `NaN` clamps to 0.
- JPEG is written at quality 90 — a documented constant rather than a silent one. Use PNG
  when the bytes must survive.
- An `Image3D` is declined (the expression stays unevaluated): a volume has no single
  raster, and quietly writing its middle slice would misreport what was exported.
- Writing is by the vendored `stb_image_write` (public domain).

#### Basic Examples

```mathematica
In[1]:= img = Image[Table[{N[i/24], N[j/32], N[Mod[i + j, 8]]/8}, {i, 1, 24}, {j, 1, 32}], "Real"];

In[2]:= Export["/tmp/mathilda_doc_e.png", img]
Out[2]= /tmp/mathilda_doc_e.png

In[3]:= FileExistsQ["/tmp/mathilda_doc_e.png"]
Out[3]= True

In[4]:= ImageDimensions[Import["/tmp/mathilda_doc_e.png"]]
Out[4]= {32, 24}
```

#### Scope

```mathematica
In[1]:= img = Image[Table[{N[i/24], N[j/32], N[Mod[i + j, 8]]/8}, {i, 1, 24}, {j, 1, 32}], "Real"];

In[2]:= Export["/tmp/mathilda_doc_e.jpg", img]
Out[1]= /tmp/mathilda_doc_e.jpg

In[2]:= Export["/tmp/mathilda_doc_e.bmp", img]
Out[2]= /tmp/mathilda_doc_e.bmp

In[3]:= Export["/tmp/mathilda_doc_e.tga", img]
Out[3]= /tmp/mathilda_doc_e.tga

In[4]:= (* a grey image writes a 1-channel file *)
ImageChannels[Import[Export["/tmp/mathilda_doc_eg.png", Image[Table[N[i/16], {i, 1, 16}, {j, 1, 16}], "Real"]]]]
Out[4]= 1
```

#### Options

```mathematica
In[1]:= img = Image[Table[{N[i/24], N[j/32], N[Mod[i + j, 8]]/8}, {i, 1, 24}, {j, 1, 32}], "Real"];

In[2]:= (* the format may be stated rather than inferred, which is the only way to write a file with no extension *)
Export["/tmp/mathilda_doc_noext", img, "PNG"]
Out[1]= /tmp/mathilda_doc_noext

In[2]:= ImageDimensions[Import["/tmp/mathilda_doc_noext", "Image"]]
Out[2]= {32, 24}
```

#### Properties & Relations

```mathematica
In[1]:= img = Image[Table[{N[i/24], N[j/32], N[Mod[i + j, 8]]/8}, {i, 1, 24}, {j, 1, 32}], "Real"];

In[2]:= (* out-of-range samples clamp to the ends rather than wrapping *)
ImageData[Import[Export["/tmp/mathilda_doc_clamp.png", Image[{{2.0, -1.0}, {1.0, 0.0}}, "Real"]]]]
Out[1]= {{1.0, 0.0}, {1.0, 0.0}}

In[2]:= (* a volume is declined rather than silently reduced to a slice *)
Head[Export["/tmp/mathilda_doc_vol.png", Image3D[Table[0.5, {z, 1, 2}, {y, 1, 2}, {x, 1, 2}], "Real"]]]
Out[2]= Export
```

## Put
Writes one or more expressions to a file, replacing any prior contents.
- `expr >> "filename"` — shorthand for `Put[expr, "filename"]`.
- `expr >> filename` — bare-word filename, equivalent to `expr >> "filename"`. The bare name follows the same token rules as `<<file` (see `Get`): `x >> f[1]` writes to `"f[1]"`, and `x >> f+g` is a syntax error.
- `Put[expr, "filename"]`
- `Put[expr_1, expr_2, ..., expr_n, "filename"]` — writes the expressions one per line.
- `Put["filename"]` — creates an empty file (or truncates an existing one).

**Features**:
- `Protected`.
- The last argument must be a string; it is interpreted as a filename.
- Each `expr_i` is rendered with the standard printer (the same form used at the REPL) and followed by a single `\n`.
- Truncates the file before writing — preserves nothing from a prior `Put`/`PutAppend`.
- Returns `Null` on success and `$Failed` on I/O error.

**Parser notes**:
- `>>` has a low precedence (`30`) — below `Set`/`SetDelayed` (`40`) and above `CompoundExpression` (`10`). `a + b >> "f"` therefore parses as `Put[a + b, "f"]`.
- The bare-word filename accepts identifier characters plus `.`, `/`, `\`, `-`, `_`, `~`, and `$`.

## PutAppend
Like `Put`, but appends to the file rather than truncating it.
- `expr >>> "filename"` — shorthand for `PutAppend[expr, "filename"]`.
- `expr >>> filename` — bare-word filename, equivalent to `expr >>> "filename"`.
- `PutAppend[expr, "filename"]`
- `PutAppend[expr_1, ..., expr_n, "filename"]`

**Features**:
- `Protected`.
- Creates the file if it does not exist; otherwise preserves prior contents and appends new lines.
- Returns `Null` on success and `$Failed` on I/O error.

**Example**:
```
FactorInteger[40320]      >>  "factorizations"
FactorInteger[479001600]  >>> "factorizations"
Get["factorizations"]     (* {{2, 10}, {3, 5}, {5, 2}, {7, 1}, {11, 1}} *)
```

## FileExistsQ
Tests for the existence of a filesystem object at the given path.
- `FileExistsQ["name"]` — `True` if anything (file, directory, symlink, FIFO, socket, device, ...) exists at `"name"`, `False` otherwise. `File["name"]` is accepted, and a leading `~` is the home directory.

**Features**:
- `Protected`.
- `"name"` is interpreted relative to the current working directory. `$Path` is not searched.
- Implemented with `lstat()`, so dangling symlinks count as existing.
- Leaves the call unevaluated when given the wrong arity, a symbolic argument, or any non-string atom.

## FileSize
Gives the size of a file, in bytes, as an integer.
- `FileSize["name"]` — number of bytes in the file `"name"`, as an `Integer` (not a `Quantity`).

**Features**:
- `Protected`.
- `"name"` is interpreted relative to the current working directory. `$Path` is not searched.
- Implemented with `stat()`, so symbolic links are followed and the target file's size is reported.
- Returns `$Failed` and prints a `FileSize::nffil` message when the file cannot be found. The message respects `Quiet[]` and is visible to `Check[]`.
- Leaves the call unevaluated when given the wrong arity, a symbolic argument, or any non-string atom.

#### Basic Examples

```mathematica
In[1]:= Put[2, 3, 5, 7, 11, "/tmp/mathilda_io_fs.txt"]  (* writes "2\n3\n5\n7\n11\n" = 11 bytes *)

In[2]:= FileSize["/tmp/mathilda_io_fs.txt"]
Out[2]= 11

In[3]:= Head[FileSize["/tmp/mathilda_io_fs.txt"]]  (* a plain Integer, not a Quantity *)
Out[3]= Integer

In[4]:= Quiet[FileSize["/tmp/mathilda_no_such_file_xyz"]]  (* missing file -> $Failed *)
Out[4]= $Failed
```

## FileExtension
Returns the trailing file extension of a path's leaf component.
- `FileExtension["name"]` — substring after the last `.` in the leaf, with no leading dot.

**Features**:
- `Protected`.
- Pure string operation — does not touch the filesystem.
- Returns `""` when the leaf has no extension, when it ends with `.`, when the leaf has only a leading `.` (e.g. `".bashrc"`), or when the path has the form of a directory (ends with `/`).
- Always ignores everything up to and including the final `/`.

## FileBaseName
Returns the leaf component of a path with its trailing extension removed.
- `FileBaseName["name"]` — leaf minus the suffix that `FileExtension` would return (and minus the separating `.`).

**Features**:
- `Protected`.
- Pure string operation — does not touch the filesystem.
- Drops everything up to and including the final `/`.
- Only the last extension is split off: `FileBaseName["file.tar.gz"]` is `"file.tar"`.
- When the leaf has no extension, the leaf is returned verbatim (including trailing `.` or leading `.`).

**Example**:
```
FileExtension["report.tar.gz"]   (* "gz" *)
FileBaseName["report.tar.gz"]    (* "report.tar" *)
FileExtension["/etc/.bashrc"]    (* "" *)
FileBaseName["/etc/.bashrc"]     (* ".bashrc" *)
```

## DirectoryName
Returns the directory portion of a path, including its final separator.
- `DirectoryName["name"]`: everything up to and including the last separator.
- `DirectoryName["name", n]`: goes up `n` levels (`n` a positive machine integer).

**Features**:
- `Protected`.
- Pure string operation; it does not touch the filesystem.
- Runs of separators are collapsed first, and a trailing separator names a
  directory: `DirectoryName["/a/b/"]` is `"/a/"`.
- `""` when there is no directory portion (`"c.txt"`, `"/"`, or more levels than the path has).
- `DirectoryName::string` for a non-string, `DirectoryName::intpm` for a bad `n`.
  Both leave the call unevaluated.
- `\` is also a separator on Windows.

**Example**:
```
DirectoryName["/a/b/c.txt"]      (* "/a/b/" *)
DirectoryName["/a/b/c.txt", 2]   (* "/a/" *)
DirectoryName["c.txt"]           (* "" *)
DirectoryName[$InputFileName]    (* the running script's directory *)
```

## ParentDirectory
The parent of an existing directory, as an absolute canonical path.
- `ParentDirectory[]`: the parent of the working directory.
- `ParentDirectory["dir"]`: the parent of `dir`. A relative `dir` is resolved
  against the working directory, and symbolic links and `..` are resolved.
- `ParentDirectory["dir", n]`: goes `n` levels up. The root is its own parent.

**Features**:
- `Protected`.
- A path that is not an existing directory emits `ParentDirectory::dirnex` and
  leaves the call unevaluated. This differs from `DirectoryName`, which is a
  pure string operation.

## FileNameTake
Selects components of a file name, as `FileNameSplit` gives them, and rejoins them.

| Call | Result for `"/a/b/c.txt"` |
|---|---|
| `FileNameTake[name]` | `"c.txt"` (the last component) |
| `FileNameTake[name, n]` / `-n` | first / last `n`: `2` → `"/a"`, `-2` → `"b/c.txt"` |
| `FileNameTake[name, {k}]` | the `k`-th: `{2}` → `"a"`, `{1}` → `"/"` |
| `FileNameTake[name, {m, n}]` | `{2, 3}` → `"a/b"` |

- An `n` larger than the number of components takes them all.
- `0`, or an index `{k}` out of range, gives `""`.
- Any other spec emits `FileNameTake::seqso`.
- `Protected`.

## AbsoluteFileName
`AbsoluteFileName["name"]` gives the absolute canonical path of an existing file
or directory. A missing one emits `AbsoluteFileName::fdnfnd` and gives
`$Failed`. `Protected`.

## FileType
`FileType["name"]` gives `File`, `Directory`, or `None` if nothing exists at
`name`. Symbolic links are followed; `File["name"]` and a leading `~` are
accepted. Any other argument emits `FileType::badfile`. `Protected`.

## DirectoryQ
`DirectoryQ["name"]` gives `True` if `name` is an existing directory (links
followed; `File["name"]` and a leading `~` accepted), and `False` otherwise --
including for a file, `""`, a list, or any non-string argument. A wrong
argument count emits `DirectoryQ::argx` and leaves the call unevaluated.
`Protected`.

## ToFileName
The legacy file-name builder, superseded by `FileNameJoin`.
- `ToFileName[{"dir1", ...}, "name"]` gives `"dir1/.../name"`.
- `ToFileName[dirs]` gives the directory with a trailing separator:
  `ToFileName[{"a", "b"}]` is `"a/b/"`.
- `dirs` may be a single string.
- Doubled separators are collapsed (`{"a/", "b"}` gives `"a/b/"`), and empty
  pieces are skipped.
- `Protected`.

## DeleteFile
`DeleteFile["file"]` or `DeleteFile[{"file1", ...}]` deletes files and returns `Null`.
- A missing file emits `DeleteFile::fdnfnd`, and a directory emits
  `DeleteFile::fdir`. In either case the result is `$Failed`, and the other
  files of a list are still deleted.
- `Protected`.

## FileNames
Lists files by name pattern.
- `FileNames[]`: every entry in the working directory.
- `FileNames[forms]`: entries whose names match `forms`.
- `FileNames[forms, dirs]`: searches the directory, or list of directories, `dirs`.
- `FileNames[forms, dirs, n]`: also searches subdirectories down to depth `n`
  (a non-negative integer or `Infinity`).

**Features**:
- `Protected`. Option `IgnoreCase -> Automatic`, which is case-insensitive on
  macOS and Windows and case-sensitive elsewhere. `True` or `False` overrides it.
- `forms` is a string pattern (`"*.m"`, with `*` and `@`), a `StringExpression`,
  a `RegularExpression`, or a list of these. It is matched against each
  entry's name.
- A string form containing a separator moves its directory part onto the
  search directory: `FileNames["sub/*.m"]` is `{"sub/d.m"}`.
- With `dirs` given, each name carries its directory (`"sub/d.m"`, `"./a.m"`).
  Without `dirs`, names are bare.
- Hidden files and directories count as entries. Subdirectories are found with
  `lstat`, so a symlinked directory is listed but not followed.
- The result is canonically sorted with duplicates removed. A missing
  directory contributes nothing.
- Depth `0` follows WL: `dirs` is ignored and the working directory is listed.
  A negative depth emits `FileNames::innf` and leaves the call unevaluated.
- Needs the PCRE2 string-pattern engine, like `StringMatchQ`.

**Example**:
```
FileNames["*.m", "src/internal"]
FileNames["*.m", ".", Infinity]       (* recursive *)
FileNames[{"*.c", "*.h"}, "src"]
```

## Directory / SetDirectory / ResetDirectory / DirectoryStack
The working directory and a stack of the ones left behind.
- `Directory[]`: the current working directory, as an absolute path.
- `SetDirectory["dir"]`: makes `dir` the working directory, pushes the previous
  one onto the stack, and returns the new directory in canonical form
  (`SetDirectory["./a/../a/"]` gives `".../a"`).
- `SetDirectory[]`: goes to `$HomeDirectory`.
- `ResetDirectory[]`: pops the stack and returns to that directory.
- `DirectoryStack[]`: the stack, most recent first. Arguments are ignored, as in WL.

**Features**:
- All `Protected`.
- Relative paths everywhere (`FileNames`, `Get`, `OpenRead`, ...) resolve
  against the current directory.
- If `dir` is missing or is not a directory, `SetDirectory::cdir` is emitted and
  the result is `$Failed`. A non-string argument emits `SetDirectory::badfile`
  and leaves the call unevaluated. In both cases the stack is unchanged.
- `ResetDirectory[]` on an empty stack emits `ResetDirectory::dtop` and returns
  the current directory.
- `getcwd`/`chdir` (`_getcwd`/`_chdir` on Windows).

```
SetDirectory["src"]; FileNames["*.c"]; ResetDirectory[]
```

## $HomeDirectory
Your "home" directory, as the full name of the directory as a string; on
multi-user systems, the main directory of the current user. It is read from
`HOME` at startup; on Windows from `USERPROFILE`, or else `HOMEDRIVE` +
`HOMEPATH`. `Protected`.

```
$HomeDirectory   (* "/Users/user" *)
```

## $InstallationDirectory
The top-level directory in which the Mathilda installation resides, as the full
(absolute, symlink-resolved) name of the directory as a string. Implemented in
`src/loadmodule.c` (`mathilda_installation_directory`) and fixed at startup.
It is read off the directory that actually serves the bootstrap `init.m`, so it
always names the tree the running process uses, whatever the current directory
or `$MATHILDA_HOME`:
- `<root>/src/internal/` → `<root>` (a binary run from its source tree);
- `<prefix>/share/mathilda/internal/` → `<prefix>` (an installed build);
- `<dir>/internal/` → `<dir>`; with nothing resolvable, the executable's directory.

`Protected`.

```
$InstallationDirectory   (* "/Users/user/Mathilda_dev/Mathilda" *)
```

## Environment
The value of an operating-system environment variable.
- `Environment["var"]`: the value of `var` as a string, or `$Failed` if it is not set.

**Features**:
- `Protected`.
- A non-string argument emits `Environment::string` and leaves the call
  unevaluated. A wrong argument count emits `Environment::argx`.

```
Environment["HOME"]          (* "/Users/..." *)
Environment["NO_SUCH_VAR"]   (* $Failed *)
```

## $TemporaryDirectory
The directory for temporary files: `TMPDIR`, else (Windows) `TEMP` or `TMP`,
else `/tmp`. A trailing separator is dropped, so
`FileNameJoin[{$TemporaryDirectory, "f"}]` never doubles it. `Protected`.

## $BaseDirectory
The base directory in which system-wide files to be loaded by Mathilda are
conventionally placed (configuration files, packages, kernel init files), as a
string. The directory need not exist; it is a convention, not a guarantee.

| Platform | `$BaseDirectory` |
|----------|------------------|
| macOS | `/Library/Mathilda` |
| Linux / Unix | `/usr/share/Mathilda` |
| Windows | `%PROGRAMDATA%\Mathilda` (`C:\ProgramData\Mathilda` when unset) |

The `MATHILDA_BASE` environment variable, when set, replaces it (a trailing
separator is dropped). `Protected`.

## $UserBaseDirectory
The per-user counterpart of `$BaseDirectory`, as a string.

| Platform | `$UserBaseDirectory` |
|----------|----------------------|
| macOS | `~/Library/Mathilda` |
| Linux / Unix | `~/.Mathilda` |
| Windows | `%APPDATA%\Mathilda` (else `<home>\AppData\Roaming\Mathilda`) |

`~` is `$HomeDirectory`. The `MATHILDA_USERBASE` environment variable, when set,
replaces it. Left unset when there is no home directory. `Protected`.

```
{$BaseDirectory, $UserBaseDirectory}
(* {"/Library/Mathilda", "/Users/user/Library/Mathilda"} *)
```

## $PreferencesDirectory
The directory in which user preferences are stored, as a string. As on the
reference system on every platform, it is the per-user base directory, so it
equals `$UserBaseDirectory` (`~/Library/Mathilda` on macOS, `~/.Mathilda` on
Unix, `%APPDATA%\Mathilda` on Windows) and follows `MATHILDA_USERBASE`.
`Protected`.

## $UserDocumentsDirectory
The user's default directory for documents, as a string:
- `XDG_DOCUMENTS_DIR`, when set (Linux / Unix only, the freedesktop user-dirs
  convention);
- else `~/Documents` (`<home>\Documents` on Windows), when that directory exists;
- else `$HomeDirectory`.

Left unset when there is no home directory. `Protected`.

```
{$PreferencesDirectory, $UserDocumentsDirectory}
(* {"/Users/user/Library/Mathilda", "/Users/user/Documents"} *)
```

## $PathnameSeparator
The separator between directory names in a file name: `"/"`, or `"\\"` on
Windows. It is the default separator `FileNameJoin` uses. `Protected`.

## $InputFileName
The absolute path of the file currently being read, or `""` when none is.
- Set while a `-file` script runs, and while `Get` (or an internal module load)
  reads a file. It is restored afterwards, so a `Get` inside a script sees the
  inner file, and the script sees its own path again once the `Get` returns.
- The path is made absolute with `realpath` (`_fullpath` on Windows). If that
  fails, the path is used as given.
- No attributes, as in WL.

## FileNameJoin
Assembles a file name from a list of path components (or canonicalizes a lone name).
- `FileNameJoin[{"n1", "n2", ...}]` — join the components with the OS pathname separator.
- `FileNameJoin["name"]` — canonicalize a single name, normalizing separators to the OS form.
- `FileNameJoin[..., OperatingSystem->"os"]` — select the separator convention; `"os"` is `"Windows"`, `"MacOSX"`, or `"Unix"`.

**Features**:
- `Protected`.
- Pure string operation — does not touch the filesystem.
- Components may themselves contain separators; each is split into segments and rejoined, so duplicate and trailing separators collapse (`{"a//b", "c"}` → `"a/b/c"`).
- An empty (or separator-led) leading component yields an absolute path: `{"", "usr", "bin"}` → `"/usr/bin"`.
- `"Windows"` uses `\` and preserves a leading `\\server\share` UNC prefix as a single unit; `"MacOSX"`/`"Unix"` use `/`. The default is the host operating system's separator.
- `Options[FileNameJoin]` reports the `OperatingSystem` default.
- `FileNameJoin[]` prints `FileNameJoin::argx` and stays unevaluated; a non-string/non-list argument, a list containing a non-string, or an unknown OS leaves the call unevaluated.

#### Basic Examples

```mathematica
In[1]:= FileNameJoin[{"usr", "local", "bin"}]
Out[1]= "usr/local/bin"

In[2]:= FileNameJoin[{"a/b", "c"}]  (* components may themselves contain separators *)
Out[2]= "a/b/c"

In[3]:= FileNameJoin[{"", "etc", "hosts"}]  (* an empty leading part gives an absolute path *)
Out[3]= "/etc/hosts"

In[4]:= FileNameJoin[{"dir1", "dir2"}, OperatingSystem -> "Windows"]  (* backslash separator *)
Out[4]= "dir1\dir2"
```

## FileNameSplit
Splits a file name into the `List` of its path components — the structural inverse of `FileNameJoin`.
- `FileNameSplit["name"]` — split using the OS pathname separator.
- `FileNameSplit[..., OperatingSystem->"os"]` — select the separator convention; `"os"` is `"Windows"`, `"MacOSX"`, or `"Unix"`.

**Features**:
- `Protected`.
- Pure string operation — does not touch the filesystem.
- A leading pathname separator marks an absolute path and yields a leading `""` part; trailing and duplicate separators are dropped (`"a//b/"` → `{"a", "b"}`).
- `"Windows"` treats a leading `\\server\share` UNC prefix as a single part and a drive like `C:` as an ordinary first part; `"MacOSX"`/`"Unix"` split on `/`. The default is the host operating system's separator.
- `Options[FileNameSplit]` reports the `OperatingSystem` default.
- `FileNameJoin[FileNameSplit[name]]` reconstructs a canonicalized `name`.
- `FileNameSplit[]` prints `FileNameSplit::argx` and stays unevaluated; a non-string argument or an unknown OS leaves the call unevaluated.

#### Basic Examples

```mathematica
In[1]:= FileNameSplit["/home/user/data.csv"]  (* leading "" marks an absolute path *)
Out[1]= {"", "home", "user", "data.csv"}

In[2]:= FileNameSplit["a//b/c/"]  (* duplicate and trailing separators are dropped *)
Out[2]= {"a", "b", "c"}

In[3]:= FileNameJoin[FileNameSplit["/usr/local/bin"]]  (* Join inverts Split *)
Out[3]= "/usr/local/bin"

In[4]:= FileNameSplit["C:\\path\\file.txt", OperatingSystem -> "Windows"]  (* a drive is an ordinary first part *)
Out[4]= {"C:", "path", "file.txt"}
```

## FilePrint
Prints the raw textual contents of a file to standard output.
- `FilePrint["file"]` — print every line.
- `FilePrint["file", n]` (n > 0) — print the first `n` lines.
- `FilePrint["file", -n]` (n > 0) — print the last `n` lines.
- `FilePrint["file", m;;n]` — print lines `m` through `n` inclusive.
- `FilePrint["file", m;;n;;s]` — same, with step `s` (positive or negative).

**Features**:
- `Protected`.
- Bytes pass through verbatim via `fwrite`, including embedded NULs and non-UTF-8 sequences.
- Lines are 1-indexed; negative indices inside the `Span` count from the end (`-1` is the last line).
- `All` may appear in any `Span` slot (`All;;-1`, `1;;All;;2`, ...) and resolves to that slot's natural endpoint.
- A positive integer larger than the file's line count clamps to "print everything"; the same applies to `-n`.
- When the file's final line lacks a trailing `\n` and the selection actually emits it, `FilePrint` adds one so the next REPL prompt isn't appended to the file content.
- Bad selectors (zero step, wrong types, wrong arity) leave the call unevaluated rather than producing partial output.
- Returns `Null` on success and `$Failed` (with a `FilePrint::noopen` diagnostic) when the file cannot be opened.

**Example**:
```
FilePrint["/etc/hosts"]                  (* whole file *)
FilePrint["/etc/hosts", 3]               (* first 3 lines *)
FilePrint["/etc/hosts", -3]              (* last 3 lines  *)
FilePrint["/etc/hosts", 2;;5]            (* lines 2..5    *)
FilePrint["/etc/hosts", 5;;1;;-1]        (* lines 5..1 reversed *)
FilePrint["/etc/hosts", 1;;-1;;2]        (* every other line *)
```
