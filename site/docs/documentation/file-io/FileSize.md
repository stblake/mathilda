# FileSize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FileSize["name"]`**

gives the number of bytes in the file with the specified name.

<details>
<summary>Notes</summary>

In FileSize\["name"\], name is interpreted relative to your current directory. FileSize does not search $Path. FileSize follows symbolic links, reporting the size of the target file. FileSize gives the size as an integer count of bytes, not a Quantity. FileSize returns $Failed and prints a message if the file cannot be found.

</details>

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`.
- `"name"` is interpreted relative to the current working directory. `$Path` is not searched.
- Implemented with `stat()`, so symbolic links are followed and the target file's size is reported.
- Returns `$Failed` and prints a `FileSize::nffil` message when the file cannot be found. The message respects `Quiet[]` and is visible to `Check[]`.
- Leaves the call unevaluated when given the wrong arity, a symbolic argument, or any non-string atom.

**Attributes:** `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_files.c`](https://github.com/stblake/mathilda/blob/main/tests/test_files.c)
