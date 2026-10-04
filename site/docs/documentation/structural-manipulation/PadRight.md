# PadRight

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PadRight[list, n]`**

makes a list of length n by padding list with zeros on the right.

**`PadRight[list, n, x]`**

pads by repeating the element x.

**`PadRight[list, n, {x1, x2, ...}]`**

pads by cyclically repeating the elements xi.

**`PadRight[list, n, padding, m]`**

leaves a margin of m elements of padding on the left.

**`PadRight[list, {n1, n2, ...}]`**

makes a nested list with length ni at level i.

**`PadRight[list]`**

pads a ragged array list with zeros to make it full.

<details>
<summary>Notes</summary>

A negative length pads on the left; a negative margin truncates leading elements. The head of list need not be List.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= PadLeft[{a, b, c}, 10, {x, y, z}]
Out[1]= {z, x, y, z, x, y, z, a, b, c}

In[2]:= PadRight[{{a, b}, {c}}, {3, 5}]
Out[2]= {{a, b, 0, 0, 0}, {c, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}
```

### Applications (4)

Pad with zeros on the right to length 5

```mathematica
In[3]:= PadRight[{a, b, c}, 5]
Out[3]= {a, b, c, 0, 0}
```

A dimension list builds a full rectangular array

```mathematica
In[4]:= PadRight[{{a, b}, {c}}, {3, 5}]
Out[4]= {{a, b, 0, 0, 0}, {c, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}
```

Automatic pads a ragged array to rectangular

```mathematica
In[5]:= PadRight[{{a, b}, {c}}]
Out[5]= {{a, b}, {c, 0}}
```

A shorter length keeps the leading elements

```mathematica
In[6]:= PadRight[{a, b, c}, 2]
Out[6]= {a, b}
```

## Implementation notes

**Algorithm.** `PadRight` is the exact mirror of `PadLeft`, sharing the recursive
`pr_build` engine with `pad_left == 0`. `pad_dispatch` parses the same specs: an
`Integer` length `n`, a `List {n1, ..., nk}` of per-level lengths (building a full
nested array, with a nested padding block tiled), or `Automatic`/omitted (pad a
ragged array to full rectangular, dimensions found by `pr_scan_dims`). The
default padding is the `Integer` `0`; a single element repeats, a `List` repeats
cyclically (via `pr_pad_at`'s floored-modulo indexing); the margin argument
leaves padding on the *left* (`PadRight`) and a negative margin truncates leading
elements. For a non-negative length `PadRight` places the original row
left-aligned and pads on the right; a negative length pads on the left. The head
of `list` need not be `List`.

**Data structures.** Identical to `PadLeft`: a heap `int64_t* dimv` of per-level
lengths, an `int64_t* coords` path threaded down the recursion to drive the
cyclic padding, and per-level source rows read directly from the input args;
`pr_is_atomic` keeps the scan and builder out of `Rational`/`Complex` nodes.

**Complexity / limits.** `O(total output elements)`. A rank-1 visible/packed
`NDArray` with a compatible fill (or pure truncation) takes the native buffer
path `ndstruct_pad` — `PadRight` is on `pack.c`'s `AWARE` and `INT64_OK` lists;
otherwise the array is materialised once and the `List` engine runs (it does not
repack, because padding an exact `0` into a `float64` buffer gives a mixed list).
There is **no** `Compile[]` lowering (`CompileDiagnostics` reports
`Compiled -> False`). `Protected`.

**Attributes:** `Protected`.

## References

**See also:** [PadLeft](../../structural-manipulation/PadLeft/), [List](../../other-advanced/List/)

- Source: [`src/list/pad.c`](https://github.com/stblake/mathilda/blob/main/src/list/pad.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_padright.c`](https://github.com/stblake/mathilda/blob/main/tests/test_padright.c)

## Notes & additional examples

### Notes

`PadRight` is the exact mirror of `PadLeft`, padding on the right instead of the
left. `PadRight[list, n]` pads with `0`; `PadRight[list, n, x]` repeats `x`, and a
list padding is tiled cyclically. A length shorter than the list keeps its
*first* `n` elements; a negative length pads on the left. `PadRight[list, {n1,
n2, ...}]` builds a full nested array, and `PadRight[list]` (or `Automatic`) pads
a ragged array to the smallest enclosing rectangle. The head of `list` need not
be `List`. A rank-1 packed/`NDArray` buffer takes a native fast path; padding an
exact `0` into a float buffer yields a mixed list that cannot repack. There is no
`Compile[]` lowering.
