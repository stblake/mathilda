# FilePrint

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FilePrint["file"]
    prints out the raw textual contents of file.
FilePrint["file", n]
    prints out the first n raw textual lines of file.
FilePrint["file", -n]
    prints out the last n raw textual lines of file.
FilePrint["file", m;;n]
    prints out lines m through n of file.
FilePrint["file", m;;n;;s]
    prints out lines m through n of file in steps of s; s may be negative
    to walk the range backwards.
FilePrint returns Null on success and $Failed if the file cannot be opened.
Negative indices inside the Span count from the end of the file (-1 is the last line).
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`.
- Bytes pass through verbatim via `fwrite`, including embedded NULs and non-UTF-8 sequences.
- Lines are 1-indexed; negative indices inside the `Span` count from the end (`-1` is the last line).
- `All` may appear in any `Span` slot (`All;;-1`, `1;;All;;2`, ...) and resolves to that slot's natural endpoint.
- A positive integer larger than the file's line count clamps to "print everything"; the same applies to `-n`.
- When the file's final line lacks a trailing `\n` and the selection actually emits it, `FilePrint` adds one so the next REPL prompt isn't appended to the file content.
- Bad selectors (zero step, wrong types, wrong arity) leave the call unevaluated rather than producing partial output.
- Returns `Null` on success and `$Failed` (with a `FilePrint::noopen` diagnostic) when the file cannot be opened.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
