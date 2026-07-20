# FileNameJoin

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FileNameJoin[{"name1", "name2", ...}]
    joins the namei into a file name suitable for your current operating system.
The namei can be individual names or file paths containing pathname separators.
FileNameJoin[{"", "name1", ...}] gives an absolute file path beginning with a pathname separator.
FileNameJoin["name"] canonicalizes name, making pathname separators appropriate for your operating system.
FileNameJoin[..., OperatingSystem->"os"] yields a file name in the format for the specified operating system; possible choices are "Windows", "MacOSX", and "Unix".
FileNameJoin just assembles a file name; it does not search for the file specified.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`.
- Pure string operation — does not touch the filesystem.
- Components may themselves contain separators; each is split into segments and rejoined, so duplicate and trailing separators collapse (`{"a//b", "c"}` → `"a/b/c"`).
- An empty (or separator-led) leading component yields an absolute path: `{"", "usr", "bin"}` → `"/usr/bin"`.
- `"Windows"` uses `\` and preserves a leading `\\server\share` UNC prefix as a single unit; `"MacOSX"`/`"Unix"` use `/`. The default is the host operating system's separator.
- `Options[FileNameJoin]` reports the `OperatingSystem` default.
- `FileNameJoin[]` prints `FileNameJoin::argx` and stays unevaluated; a non-string/non-list argument, a list containing a non-string, or an unknown OS leaves the call unevaluated.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
