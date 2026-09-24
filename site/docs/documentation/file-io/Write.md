# Write

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Write[stream, expr1, expr2, ...]`**

writes the expressions to an output stream in input form, followed by a newline. The stream may be an OutputStream, a "file", or File\["file"\]; a named file that is not already open is opened for writing and left open.

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`. Return `$Failed` if the file cannot be opened.
- `Write` evaluates its expression arguments before writing (use `Hold[...]` to write an unevaluated form), and output is flushed after each call so it round-trips with `Read`/`ReadList`.

**Attributes:** `Protected`.

## References

**See also:** [WriteString](../../file-io/WriteString/), [OutputStream](../../other-advanced/OutputStream/), [Read](../../file-io/Read/), [ReadList](../../file-io/ReadList/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)
