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

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= ByteCount[5]
Out[1]= 48

In[2]:= ByteCount[{1, 2, 3}]
Out[2]= 269

In[3]:= ByteCount[x^2 + 1]
Out[3]= 381
```

### Notes

`ByteCount` reports the number of bytes Mathilda uses internally to store an expression, including every subexpression node. The exact figures are an implementation detail of the build and are useful mainly for comparing the relative size of expressions.
