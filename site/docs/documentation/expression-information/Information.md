# Information

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Information[symbol] or ?symbol returns information on symbol.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_information` (`src/core.c`) looks up the symbol's docstring with `symtab_get_docstring` and returns it as a string. If none exists it returns a string `No information available for symbol "..."` using `context_display_name` for the shortened name. (The interactive `?name` syntax routes to the same docstring store.)

**Attributes:** `HoldAll`, `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_blas.c`](https://github.com/stblake/mathilda/blob/main/tests/test_blas.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_eigen.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eigen.c)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)
