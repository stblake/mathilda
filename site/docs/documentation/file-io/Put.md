# Put

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Put[expr, "filename"] or expr >> "filename"
    writes expr to the file, replacing any prior contents.
Put[expr1, expr2, ..., "filename"]
    writes a sequence of expressions to the file, each followed by a newline.
Put["filename"]
    creates an empty file with the specified name (or truncates an existing one).
expr >> "filename" is equivalent to Put[expr, "filename"]; the bare-word form
expr >> filename is equivalent to expr >> "filename".
Returns Null on success and $Failed if the file cannot be opened.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`.
- The last argument must be a string; it is interpreted as a filename.
- Each `expr_i` is rendered with the standard printer (the same form used at the REPL) and followed by a single `\n`.
- Truncates the file before writing — preserves nothing from a prior `Put`/`PutAppend`.
- Returns `Null` on success and `$Failed` on I/O error.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
