# ByteCount

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ByteCount[expr] gives the number of bytes used internally by Mathilda to store expr.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= ByteCount[5]
Out[1]= 48

In[2]:= ByteCount[{1, 2, 3}]
Out[2]= 269

In[3]:= ByteCount[x^2 + 1]
Out[3]= 381

In[4]:= ByteCount[Range[100]]
Out[4]= 5701

In[5]:= ByteCount[Factorial[50]]
Out[5]= 48
```

## Implementation notes

`builtin_bytecount` (`src/core.c`) recursively sums `byte_count_internal` over the tree: `sizeof(Expr)` per node, plus `strlen+1` for symbol/string payloads and `sizeof(Expr*) * arg_count` for each function's argument array, descending into the head and all arguments. It returns an integer; the count is a structural estimate and does not account for GMP/MPFR limb storage of bigints/reals.

- `Protected`.
- Uses `sizeof()` in C and measures the internal AST memory allocation boundaries, dynamically capturing sizes of individual strings, symbols, allocated blocks, arrays, and expression structs.
- Counts the payload of leaf atoms that own out-of-node storage: `EXPR_BIGINT` (GMP limbs), `EXPR_NDARRAY` (the `dims[]` array plus the flat data buffer, sized by element count and dtype width), and `EXPR_MPFR` (significand storage, scaling with precision). For an `NDArray`, the buffer dominates, so `ByteCount` scales with the number of elements and the dtype's bytes-per-element.

**Attributes:** `Protected`.

## References

**See also:** [NDArray](../../linear-algebra/NDArray/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)

## Notes & additional examples

### Notes

`ByteCount` reports the number of bytes Mathilda uses internally to store an expression, including every subexpression node. The exact figures are an implementation detail of the build and are useful mainly for comparing the relative size of expressions. Note that an arbitrary-precision integer such as `Factorial[50]` reports the same node size as a small integer (`48`) because the count measures the expression node, not the heap-allocated GMP limbs behind a bignum — whereas a 100-element list scales linearly with its elements.
