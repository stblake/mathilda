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

**Algorithm.** `builtin_fileprint` slurps the whole file into one buffer (`fileprint_slurp`, binary mode via `fopen`/`ftell`/`fread`); on open failure it prints `FilePrint::noopen` and returns `$Failed`. `fileprint_split_lines` then records one `LineSpan {start, len}` per line in a single pass (each line's trailing `\n` is included in its length, and an unterminated final line is captured). The optional selector is decoded by `fileprint_decode_selector` into a 1-based `(start, end, step)` triple: a positive integer `n` selects the first `n` lines, negative the last `|n|`, and `Span[m, n]` / `Span[m, n, s]` an explicit (possibly reverse, signed-step) range with `All` and negative endpoints resolved by `fileprint_resolve_span_slot`. `fileprint_emit` writes the selected lines verbatim with `fwrite` (preserving embedded NULs / non-UTF-8 bytes), appending a `\n` only if the last emitted line lacked one. Returns the symbol `Null`; an ill-formed selector leaves the call unevaluated rather than printing partially. `ATTR_PROTECTED`.

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

- Source: [`src/files.c`](https://github.com/stblake/mathilda/blob/main/src/files.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Put[x^2 + 1, "/tmp/mathilda_demo.m"]
Out[1]= Null

In[2]:= FilePrint["/tmp/mathilda_demo.m"]
1 + x^2
Out[2]= Null
```

### Notes

`FilePrint["file"]` writes the raw textual contents of the file to stdout and returns `Null`; it does not parse or evaluate the file (use `Get` for that).
