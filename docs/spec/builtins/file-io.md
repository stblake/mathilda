# File I/O

Builtins implemented in `src/readwrite.c` (`Get`/`Put`/`PutAppend`) and `src/files.c` (`FileExistsQ`, `FileExtension`, `FileBaseName`).

## Get
Reads a sequence of Mathilda expressions from a file, evaluates each in order, and returns the value of the last one.
- `Get["filename"]`

**Features**:
- `Protected`.
- Returns `$Failed` if the file cannot be opened.
- Used by the REPL bootstrap to load `src/internal/init.m` (and the rules it pulls in).
- Files conventionally end with `.m`.

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
