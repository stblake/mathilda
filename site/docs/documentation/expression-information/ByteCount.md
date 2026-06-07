# ByteCount

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ByteCount[expr] gives the number of bytes used internally by Mathilda to store expr.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_bytecount` (`src/core.c`) recursively sums `byte_count_internal` over the tree: `sizeof(Expr)` per node, plus `strlen+1` for symbol/string payloads and `sizeof(Expr*) * arg_count` for each function's argument array, descending into the head and all arguments. It returns an integer; the count is a structural estimate and does not account for GMP/MPFR limb storage of bigints/reals.

- `Protected`.
- Uses `sizeof()` in C and measures the internal AST memory allocation boundaries, dynamically capturing sizes of individual strings, symbols, allocated blocks, arrays, and expression structs.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
