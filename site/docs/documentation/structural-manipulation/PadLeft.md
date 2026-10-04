# PadLeft

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PadLeft[list, n]`**

makes a list of length n by padding list with zeros on the left.

**`PadLeft[list, n, x]`**

pads by repeating the element x.

**`PadLeft[list, n, {x1, x2, ...}]`**

pads by cyclically repeating the elements xi.

**`PadLeft[list, n, padding, m]`**

leaves a margin of m elements of padding on the right.

**`PadLeft[list, {n1, n2, ...}]`**

makes a nested list with length ni at level i.

**`PadLeft[list]`**

pads a ragged array list with zeros to make it full.

<details>
<summary>Notes</summary>

A negative length pads on the right; a negative margin truncates trailing elements. The head of list need not be List.

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

Pad with zeros on the left to length 5

```mathematica
In[3]:= PadLeft[{a, b, c}, 5]
Out[3]= {0, 0, a, b, c}
```

A list padding is tiled cyclically

```mathematica
In[4]:= PadLeft[{a, b, c}, 10, {x, y, z}]
Out[4]= {z, x, y, z, x, y, z, a, b, c}
```

A shorter length keeps the trailing elements

```mathematica
In[5]:= PadLeft[{1, 2, 3}, 2]
Out[5]= {2, 3}
```

A dimension list builds a full nested array

```mathematica
In[6]:= PadLeft[{{a, b}, {c}}, {3, 4}]
Out[6]= {{0, 0, 0, 0}, {0, 0, a, b}, {0, 0, 0, c}}
```

## Implementation notes

**Algorithm.** `PadLeft` and `PadRight` are exact mirrors and share one recursive
engine `pr_build`, selected by a `pad_left` flag. `pad_dispatch` parses the
length spec (arg 2): an `Integer` `n` (a single level), a `List {n1, ..., nk}`
(target length `ni` at level `i`, building a full nested array), or
`Automatic`/omitted (pad the ragged input to full rectangular — `pr_scan_dims`
records the maximum width seen at each level). Padding (arg 3) defaults to the
`Integer` `0`, may be a single element, or a `List` tiled cyclically; the margin
(arg 4) leaves that many padding elements on the far side (negative truncates).
`pr_build` lays out each level into `N = |n|` slots, placing the original row
left-aligned (for right padding) or right-aligned (for left padding) and filling
the rest from the padding block via `pr_pad_at`, which indexes into a `List`
padding with floored modulo so the tiling phase is continuous. A negative length
pads on the opposite side; the head of `list` need not be `List` and is
preserved.

**Data structures.** A heap `int64_t* dimv` of per-level target lengths, an
`int64_t* coords` path accumulated down the recursion (used to index the cyclic
padding), and the per-level source row read straight off the input node's args.
`pr_is_atomic` stops the dimension scan and the builder from descending into
`Rational`/`Complex` nodes.

**Complexity / limits.** `O(total output elements)`. A rank-1 visible/packed
`NDArray` with a dtype-compatible fill (or pure truncation) takes the native
buffer path `ndstruct_pad` — `PadLeft`/`PadRight` are on `pack.c`'s `AWARE` and
`INT64_OK` lists; otherwise the array is materialised once and the `List` engine
runs (padding an exact `0` into a `float64` buffer yields a mixed list no uniform
buffer can hold, so it deliberately does not repack). There is **no** `Compile[]`
lowering (`CompileDiagnostics` reports `Compiled -> False` for this head).
`Protected`.

**Attributes:** `Protected`.

## References

**See also:** [PadRight](../../structural-manipulation/PadRight/), [List](../../other-advanced/List/)

- Source: [`src/list/pad.c`](https://github.com/stblake/mathilda/blob/main/src/list/pad.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_padleft.c`](https://github.com/stblake/mathilda/blob/main/tests/test_padleft.c)

## Notes & additional examples

### Notes

`PadLeft[list, n]` returns a length-`n` list padded on the left with `0`;
`PadLeft[list, n, x]` repeats a given element, and `PadLeft[list, n, {x1, ...}]`
tiles a list of pad elements cyclically so the padding phase is continuous. A
length shorter than the list keeps its *last* `n` elements (left padding, so the
front is dropped); a negative length pads on the opposite side. The dimension-list
form `PadLeft[list, {n1, n2, ...}]` builds a full nested array with length `ni` at
level `i`, and the bare `PadLeft[list]` pads a ragged array to rectangular. The
head of `list` need not be `List`. A rank-1 packed/`NDArray` buffer takes a native
fast path; there is no `Compile[]` lowering.
