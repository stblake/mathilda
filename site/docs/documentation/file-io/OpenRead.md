# OpenRead

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`OpenRead["file"]`**

opens a file for reading and returns an InputStream object; also accepts File\["file"\]. Returns $Failed if the file cannot be opened.

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`. Return `$Failed` (with an `Open*::noopen` diagnostic) if the file cannot be opened.
- The integer in the returned object is an internal handle into the stream registry; the object is inert and prints as itself.
- The registry is freed on `Close` or, for anything still open, at program exit (no leaks).

**Attributes:** `Protected`.

## References

**See also:** [OpenWrite](../../file-io/OpenWrite/), [OpenAppend](../../file-io/OpenAppend/), [Close](../../file-io/Close/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)
