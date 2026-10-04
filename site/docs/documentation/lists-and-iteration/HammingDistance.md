# HammingDistance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HammingDistance[u, v]`**

Gives the number of positions at which two equal-length strings or lists differ. Returns unevaluated when the lengths differ.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= EditDistance["GGTTT", "GGGGT"]
Out[1]= 2

In[2]:= EditDistance["kitten", "sitting"]
Out[2]= 3

In[3]:= EditDistance[{1, 2, 3}, {1, 3}]
Out[3]= 1

In[4]:= HammingDistance["GGTTT", "GGGGT"]
Out[4]= 2
```

### Applications (3)

Two positions differ

```mathematica
In[5]:= HammingDistance[{1, 2, 3, 4}, {1, 0, 3, 0}]
Out[5]= 2
```

Bit-vector distance

```mathematica
In[6]:= HammingDistance[{1, 0, 1, 1}, {1, 1, 0, 1}]
Out[6]= 2
```

Compared character by character

```mathematica
In[7]:= HammingDistance["2718281828", "3141592653"]
Out[7]= 10
```

## Implementation notes

**Algorithm.** `builtin_hamming_distance` counts how many positions differ
between two equal-length sequences. After `dist_seq_pair` confirms both arguments
are strings or both lists, `dist_seq` turns each into an element array (a string
explodes to one `Expr` per byte; a list's elements are borrowed) and a single
pass tallies the positions where `expr_eq` is false. Equal lengths are required:
on a length mismatch the call is left unevaluated (Mathematica raises `::idim`
there, so declining is the faithful behaviour).

**Encoding.** Strings are compared byte by byte, so a multi-byte UTF-8 character
spans several positions — matching the ASCII/DNA use cases and noted rather than
assumed. Comparison by `expr_eq` means lists of arbitrary expressions work, not
just numbers or characters.

**Complexity / limits.** O(n) comparisons, O(n) temporary for the string case.
`ATTR_PROTECTED`. See `EditDistance` for the insert/delete/substitute metric that
does not require equal lengths.

- `Protected`.
- Elements are compared with structural equality, so the same routine serves
  strings (character by character) and lists of arbitrary expressions:
  `EditDistance[{1, 2, 3}, {1, 3}]` is `1`.
- Strings are compared **byte by byte**, so a multi-byte UTF-8 character counts
  as several elements.
- `HammingDistance` requires equal lengths and leaves the call unevaluated
  otherwise, matching Mathematica's `::idim`.
- `EditDistance` costs `O(m n)` time and `O(min(m, n))` memory (two DP rows).

**Attributes:** `Protected`.

## References

**See also:** [EditDistance](../../lists-and-iteration/EditDistance/)

- Source: [`src/list/distance.c`](https://github.com/stblake/mathilda/blob/main/src/list/distance.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)

## Notes & additional examples

### Notes

`HammingDistance[a, b]` counts the positions at which two equal-length sequences
differ. Both arguments must be strings, or both lists; a length mismatch leaves
the call unevaluated (unlike `EditDistance`, which handles differing lengths).
Comparison is by equality, so it works on bit vectors, digit strings, or lists of
arbitrary expressions. Strings are compared byte by byte, so a multi-byte UTF-8
character spans several positions.
