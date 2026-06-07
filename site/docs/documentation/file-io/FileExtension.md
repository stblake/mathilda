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

`builtin_fileextension` is pure string manipulation (no filesystem access). It isolates the filename component (after the last `/`, via `filename_component`) and finds the extension start with `extension_offset`: the offset just past the last `.` that is neither the first nor the last character of the component (so `.bashrc` and `file.` have no extension). It returns the substring from that offset — i.e. the suffix without the dot — as an `EXPR_STRING`, or `""` when there is no extension. `ATTR_PROTECTED`.

- `Protected`.
- Pure string operation — does not touch the filesystem.
- Returns `""` when the leaf has no extension, when it ends with `.`, when the leaf has only a leading `.` (e.g. `".bashrc"`), or when the path has the form of a directory (ends with `/`).
- Always ignores everything up to and including the final `/`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/files.c`](https://github.com/stblake/mathilda/blob/main/src/files.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
