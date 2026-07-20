# Get

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Get["filename"]
    reads expressions from a file, evaluates them in order, and returns the last result.
Expressions are separated by a newline or a semicolon; a trailing ; is not required.
A line break inside (), [], {} or <||> is insignificant, so one expression may span lines.
Returns $Failed if the file cannot be opened.
It is conventional to use names ending in .m for files containing Mathilda input.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_get` reads a Mathilda source file and evaluates it expression by expression, returning the last value. It opens the file (`Get::noopen` + `$Failed` on failure), slurps the entire contents into a `malloc`'d buffer, then walks the buffer with the parser's `parse_next_expression(&ptr)` — the same Pratt parser used by the REPL — `evaluate`ing each parsed expression and keeping the last non-`NULL` result (defaulting to `Null` for an empty file). Parsing stops when `parse_next_expression` returns `NULL` at end-of-input. This is the mechanism `init.m` uses to load the internal `.m` bootstrap files. `ATTR_PROTECTED`.

- `Protected`.
- Returns `$Failed` if the file cannot be opened.
- Used by the REPL bootstrap to load `src/internal/init.m` (and the rules it pulls in).
- Files conventionally end with `.m`.

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

In[2]:= Get["/tmp/mathilda_demo.m"]
Out[2]= 1 + x^2
```

### Notes

`Get["file"]` reads the expressions in a file, evaluates them in order, and returns the last result.
