# ArrayReshape

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ArrayReshape[list, dims]`**

arranges the flattened elements of list into a rectangular array of dimensions dims, dropping extra elements or padding with 0 as needed.

**`ArrayReshape[list, dims, padding]`**

uses the given padding scheme (as in ArrayPad) when list is too short.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= ArrayReshape[{a, b, c, d, e, f}, {2, 3}]
Out[1]= {{a, b, c}, {d, e, f}}

In[2]:= ArrayReshape[{1, 2, 3, 4, 5, 6, 7}, {5, 3}, x]
Out[2]= {{1, 2, 3}, {4, 5, 6}, {7, x, x}, {x, x, x}, {x, x, x}}
```

### Applications (4)

Fill a 3x4 matrix row by row

```mathematica
In[3]:= ArrayReshape[Range[12], {3, 4}]
Out[3]= {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}}
```

Too many elements: the extras are dropped

```mathematica
In[4]:= ArrayReshape[Range[6], {2, 2}]
Out[4]= {{1, 2}, {3, 4}}
```

Too few: the tail is padded with 0

```mathematica
In[5]:= ArrayReshape[{1, 2}, {2, 3}, 0]
Out[5]= {{1, 2, 0}, {0, 0, 0}}
```

Reshaping preserves the flattened order

```mathematica
In[6]:= Flatten[ArrayReshape[Range[24], {2, 3, 4}]] === Range[24]
Out[6]= True
```

## Implementation notes

**Algorithm.** `builtin_array_reshape` arranges the flattened elements of a list
into a rectangular `dims` array. `ar_parse_dims` reads the dims spec (a
non-negative integer, or a `List` of them) and computes the element total with an
overflow guard. `ar_collect` then fully flattens the input, descending `List`
heads only, into a growing buffer of borrowed leaf pointers. The output is forced
to exactly `total` owned leaves — copied/truncated when the input has enough, or
extended by `pad_scheme_extend` (default fill `0`, or a given padding scheme)
when it has too few — and `ar_build` folds that flat sequence into the nested
rectangular shape in row-major order. So up to the shared length,
`Flatten[ArrayReshape[list, dims]] == Flatten[list]`.

**Data structures / limits.** Rank is capped at `AR_MAX_RANK` (64); the leaf
buffer starts at 16 and doubles. A packed/`NDArray` first argument takes the
`ndstruct_arrayreshape` buffer fast path (a dims-header `memcpy`, no
per-element work); `ArrayReshape` is on `pack.c`'s `AWARE` list.

**Complexity / limits.** O(total) element copies plus the flatten pass.
`ATTR_PROTECTED`. A bad dims spec (empty list, negative or non-integer entry,
rank overflow, or an int64-overflowing total) or a non-list, non-NDArray first
argument leaves the call unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [ArrayPad](../../lists-and-iteration/ArrayPad/), [NDArray](../../linear-algebra/NDArray/)

- Source: [`src/list/array_reshape.c`](https://github.com/stblake/mathilda/blob/main/src/list/array_reshape.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_array_reshape.c`](https://github.com/stblake/mathilda/blob/main/tests/test_array_reshape.c)

## Notes & additional examples

### Notes

`ArrayReshape[list, dims]` lays the fully flattened elements of `list` into a
rectangular `dims` array in row-major order. If `list` has more elements than the
shape needs, the extras are dropped; if it has fewer, the tail is filled with the
padding (default `0`, or a named scheme as a third argument). Up to the shared
length, `Flatten[ArrayReshape[list, dims]] == Flatten[list]`. Dimensions must be
non-negative integers.
