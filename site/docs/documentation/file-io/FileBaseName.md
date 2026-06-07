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

- `Protected`.
- Pure string operation — does not touch the filesystem.
- Drops everything up to and including the final `/`.
- Only the last extension is split off: `FileBaseName["file.tar.gz"]` is `"file.tar"`.
- When the leaf has no extension, the leaf is returned verbatim (including trailing `.` or leading `.`).

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
