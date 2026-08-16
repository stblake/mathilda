# File I/O

Builtins implemented in `src/readwrite.c` (`Get`/`Put`/`PutAppend`) and `src/files.c` (`FileExistsQ`, `FileExtension`, `FileBaseName`, `FileNameJoin`, `FileNameSplit`, `FilePrint`).

## Get
Reads a sequence of Mathilda expressions from a file, evaluates each in order, and returns the value of the last one.
- `Get["filename"]`

**Features**:
- `Protected`.
- Returns `$Failed` if the file cannot be opened.
- Used by the REPL bootstrap to load `src/internal/init.m` (and the rules it pulls in).
- Files conventionally end with `.m`.

## LoadModule
Loads an internal Mathilda source module, resolving its location independently of
the current working directory.
- `LoadModule["relpath"]` — `relpath` is relative to the source tree's
  `src/internal` directory (e.g. `"simp/FullSimplify.m"`).

**Features**:
- `Protected`. Returns `True` if the module was located and loaded (or had already
  been loaded), `False` otherwise.
- Resolution is independent of the current working directory and tries, in
  order: `$MATHILDA_HOME/<relpath>`; `<exe_dir>/src/internal/<relpath>` and
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

## Import
Reads a raster image file and returns an `Image`.
- `Import["file"]` — decodes PNG, JPEG, BMP, GIF, TGA, PSD, HDR or PNM by content.
- `Import["file", "Image"]` — the same, stated explicitly.

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
Writes an `Image` to a raster image file.
- `Export["file", image]` — the format comes from the file extension (PNG, JPEG, BMP, TGA).
- `Export["file", image, "PNG"]` — the format stated explicitly, for a name that does not
  carry one.

**Features**:
- `Protected`.
- Returns the file name, so `Import[Export[f, img]]` is a round trip that can be written
  as a single expression.
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
- `expr >> filename` — bare-word filename, equivalent to `expr >> "filename"`.
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
- `FileExistsQ["name"]` — `True` if anything (file, directory, symlink, FIFO, socket, device, ...) exists at `"name"`, `False` otherwise.

**Features**:
- `Protected`.
- `"name"` is interpreted relative to the current working directory. `$Path` is not searched.
- Implemented with `lstat()`, so dangling symlinks count as existing.
- Leaves the call unevaluated when given the wrong arity, a symbolic argument, or any non-string atom.

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

**Example**:
```
FileNameJoin[{"dir1", "dir2", "file"}]                       (* "dir1/dir2/file" *)
FileNameJoin[{"dir1/dir2", "file"}]                          (* "dir1/dir2/file" *)
FileNameJoin[{"", "usr", "bin"}]                             (* "/usr/bin" *)
FileNameJoin[{"dir1", "dir2"}, OperatingSystem->"Windows"]   (* "dir1\dir2" *)
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

**Example**:
```
FileNameSplit["a/b/c"]                                        (* {"a", "b", "c"} *)
FileNameSplit["/home/sb/mathilda/examples/"]                  (* {"", "home", "sb", "mathilda", "examples"} *)
FileNameSplit["C:\path\file", OperatingSystem->"Windows"]     (* {"C:", "path", "file"} *)
FileNameSplit["\\server\share\path\file", OperatingSystem->"Windows"]  (* {"\\server\share", "path", "file"} *)
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
