# WriteString

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`WriteString[stream, str1, str2, ...]`**

writes the strings to an output stream with no added quotes or newline. Non-string arguments are written in input form.

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`. Return `$Failed` if the file cannot be opened.
- `Write` evaluates its expression arguments before writing (use `Hold[...]` to write an unevaluated form), and output is flushed after each call so it round-trips with `Read`/`ReadList`.

**Attributes:** `Protected`.

## References

**See also:** [Write](../../file-io/Write/), [OutputStream](../../other-advanced/OutputStream/), [Read](../../file-io/Read/), [ReadList](../../file-io/ReadList/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)
