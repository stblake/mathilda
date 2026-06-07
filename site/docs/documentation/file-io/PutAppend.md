# PutAppend

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PutAppend[expr, "filename"] or expr >>> "filename"
    appends expr to the end of the file, creating the file if it does not exist.
PutAppend[expr1, expr2, ..., "filename"]
    appends a sequence of expressions, one per line.
PutAppend works the same as Put, except that it preserves any existing
contents of the file rather than truncating them.
expr >>> filename is equivalent to expr >>> "filename".
Returns Null on success and $Failed if the file cannot be opened.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_putappend` shares `put_common` with `Put`, differing only in fopen mode: it opens the (last-argument) filename in `"a"` mode, appending each preceding expression — serialized via `expr_to_string` (the standard printer) plus a trailing `\n` — to the end of the file, creating it if absent and never truncating an existing one. `PutAppend["file"]` writes nothing. Open failure prints `PutAppend::noopen` and returns `$Failed`; success returns `Null`. `ATTR_PROTECTED`. The infix `expr >>> "file"` lowers to `PutAppend[expr, "file"]`.

- `Protected`.
- Creates the file if it does not exist; otherwise preserves prior contents and appends new lines.
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

In[2]:= PutAppend[y, "/tmp/mathilda_demo.m"]
Out[2]= Null

In[3]:= FilePrint["/tmp/mathilda_demo.m"]
1 + x^2
y
Out[3]= Null
```

### Notes

`PutAppend[expr, "file"]` (equivalently `expr >>> "file"`) adds `expr` on a new line at the end of the file, creating the file if it does not exist.
