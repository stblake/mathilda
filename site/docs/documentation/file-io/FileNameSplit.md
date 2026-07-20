# FileNameSplit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FileNameSplit["name"]
    splits a file name into a list of parts.
FileNameSplit by default uses pathname separators and other conventions suitable for your operating system.
FileNameSplit[..., OperatingSystem->"os"] uses the conventions of the specified operating system; possible choices are "Windows", "MacOSX", and "Unix".
Absolute file names that begin with a pathname separator yield a list of parts that starts with "".
Under Windows, the drive or share name is treated as the first part of the file name.
FileNameSplit just operates on names of files; it does not actually search for the file specified.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`.
- Pure string operation — does not touch the filesystem.
- A leading pathname separator marks an absolute path and yields a leading `""` part; trailing and duplicate separators are dropped (`"a//b/"` → `{"a", "b"}`).
- `"Windows"` treats a leading `\\server\share` UNC prefix as a single part and a drive like `C:` as an ordinary first part; `"MacOSX"`/`"Unix"` split on `/`. The default is the host operating system's separator.
- `Options[FileNameSplit]` reports the `OperatingSystem` default.
- `FileNameJoin[FileNameSplit[name]]` reconstructs a canonicalized `name`.
- `FileNameSplit[]` prints `FileNameSplit::argx` and stays unevaluated; a non-string argument or an unknown OS leaves the call unevaluated.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
