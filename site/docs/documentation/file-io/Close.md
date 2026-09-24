# Close

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Close[stream]`**

closes an open input or output stream and returns its file name.

**`Close["file"] closes a stream opened for that file. Returns $Failed if`**

the stream is not open.

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [InputStream](../../other-advanced/InputStream/), [OutputStream](../../other-advanced/OutputStream/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)
