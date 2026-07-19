# File I/O

Builtins implemented in `src/readwrite.c` (`Get`/`Put`/`PutAppend`) and `src/files.c` (`FileExistsQ`, `FileExtension`, `FileBaseName`, `FileNameJoin`, `FilePrint`).

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
