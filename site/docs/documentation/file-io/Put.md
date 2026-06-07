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

`builtin_put` (and `builtin_putappend`) share `put_common`, which treats the **last** argument as the filename and writes every preceding expression to it, one per line. `Put` opens the file in `"w"` mode (truncating); each `expr_i` is serialized with `expr_to_string` (the standard printer) followed by `\n`. The single-argument form `Put["file"]` runs the loop zero times, leaving an empty file. Open failure prints `Put::noopen` and returns `$Failed`; success returns the symbol `Null`. The output is intended to be re-readable with `Get`. `ATTR_PROTECTED`. The infix `expr >> "file"` lowers to `Put[expr, "file"]`.

- `Protected`.
- The last argument must be a string; it is interpreted as a filename.
- Each `expr_i` is rendered with the standard printer (the same form used at the REPL) and followed by a single `\n`.
- Truncates the file before writing — preserves nothing from a prior `Put`/`PutAppend`.
- Returns `Null` on success and `$Failed` on I/O error.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/readwrite.c`](https://github.com/stblake/mathilda/blob/main/src/readwrite.c)
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

`Put[expr, "file"]` (equivalently `expr >> "file"`) writes `expr` to the file, replacing any prior contents. Read it back with `Get`.
