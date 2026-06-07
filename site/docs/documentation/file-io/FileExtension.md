# FileExtension

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FileExtension["file"]
    gives the file extension for a file name.
FileExtension["name.ext"] gives "ext".
FileExtension gives the extension that appears after the last . in a file name.
If there are multiple endings to a file name, separated by ., FileExtension gives only the last one.
FileExtension gives "" if there is no file extension, if the file name has the form of a directory name, or ends with a . character.
FileExtension ignores any directory specification.
FileExtension by default assumes pathname separators and other conventions suitable for your operating system.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`.
- Pure string operation — does not touch the filesystem.
- Returns `""` when the leaf has no extension, when it ends with `.`, when the leaf has only a leading `.` (e.g. `".bashrc"`), or when the path has the form of a directory (ends with `/`).
- Always ignores everything up to and including the final `/`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
