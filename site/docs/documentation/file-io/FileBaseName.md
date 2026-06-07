# FileBaseName

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FileBaseName["file"]
    gives the base name for a file without its extension.
FileBaseName["name.ext"] gives "name".
FileBaseName drops all directory specifications.
FileBaseName["file.tar.gz"] gives "file.tar" — only the final extension is split off.
FileBaseName by default assumes pathname separators and other conventions suitable for your operating system.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_filebasename` is pure string manipulation (no filesystem access). It takes the filename component (everything after the last `/`, via `filename_component`), then strips a trailing extension using `extension_offset` — the suffix after the last `.` that is neither the first nor last character of the component. The base name is the component up to (not including) that `.`; with no qualifying extension the whole component (including any trailing `.`) is kept. The result is `memcpy`'d into a fresh buffer and returned as an `EXPR_STRING`. Non-string input is unevaluated. `ATTR_PROTECTED`.

- `Protected`.
- Pure string operation — does not touch the filesystem.
- Drops everything up to and including the final `/`.
- Only the last extension is split off: `FileBaseName["file.tar.gz"]` is `"file.tar"`.
- When the leaf has no extension, the leaf is returned verbatim (including trailing `.` or leading `.`).

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/files.c`](https://github.com/stblake/mathilda/blob/main/src/files.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
