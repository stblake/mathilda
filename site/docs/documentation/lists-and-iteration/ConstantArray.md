# ConstantArray

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ConstantArray[c, n]`**

generates a list of n copies of the element c.

**`ConstantArray[c, {n1, n2, ...}]`**

generates an n1 x n2 x ... nested array of copies of c.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= ConstantArray[c, 5]
Out[1]= {c, c, c, c, c}

In[2]:= ConstantArray[0, {2, 3}]
Out[2]= {{0, 0, 0}, {0, 0, 0}}
```

### Applications (3)

Any expression may be the repeated element

```mathematica
In[3]:= ConstantArray[Pi, 4]
Out[3]= {Pi, Pi, Pi, Pi}
```

A nested rectangular array

```mathematica
In[4]:= ConstantArray[1, {2, 2}]
Out[4]= {{1, 1}, {1, 1}}
```

A symbolic element is copied verbatim

```mathematica
In[5]:= ConstantArray[x, {2, 3}]
Out[5]= {{x, x, x}, {x, x, x}}
```

## Implementation notes

**Algorithm.** `builtin_constant_array` builds a flat list of `n` copies of `c`
(`ConstantArray[c, n]`), or an `n1 x ... x nk` nested array
(`ConstantArray[c, {n1, ..., nk}]`). It is `Array[]` minus the index
computation: `ca_helper` recurses to the deepest level and returns a fresh
`expr_copy(c)` at each leaf, with no indexed function call built or evaluated.
Every dimension must be a non-negative machine integer (else the call is left
unevaluated); a `0` dimension yields an empty `List` at that level. The optional
third (padding) argument is accepted for arity but not implemented, so that form
is declined.

**Buffer fast path.** When `c` is a machine number (`EXPR_INTEGER` or
`EXPR_REAL`) over a rectangular shape, the whole result is known before anything
is built, so `ndbuild_open` opens a packed `NDArray` (`NDT_INT64` or
`NDT_FLOAT64`) and the single value is written straight into the buffer — no
per-element `Expr` is allocated. `UnitVector[10^6, ...]`-scale construction thus
costs a buffer fill rather than 10⁶ boxed nodes. `ndbuild_open` declines a zero
dimension or sub-threshold size, routing those shapes back to `ca_helper`.

**Data structures / limits.** `Expr**` children for the boxed path; a dense
`int64`/`double` buffer for the packed path (rank under `NDARRAY_MAX_RANK`). A
symbolic or compound `c` always takes the boxed path. `ATTR_PROTECTED`.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [NDArrayQ](../../other-advanced/NDArrayQ/)

- Source: [`src/list/constant_array.c`](https://github.com/stblake/mathilda/blob/main/src/list/constant_array.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_characteristicpolynomial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_characteristicpolynomial.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_constant_array.c`](https://github.com/stblake/mathilda/blob/main/tests/test_constant_array.c)

## Notes & additional examples

### Notes

`ConstantArray[c, n]` is a flat list of `n` copies of `c`; `ConstantArray[c, {n1,
..., nk}]` is an `n1 x ... x nk` nested array. The element `c` is copied
verbatim, so it may be a symbol, a number, or a compound structure such as a
matrix. A dimension of `0` yields an empty list at that level; dimensions must be
non-negative machine integers.

When `c` is a machine number over a rectangular shape, the whole result is a
packed array written straight into a buffer, so building a large constant array
costs a memory fill rather than one boxed element per cell.
