# EditDistance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EditDistance[u, v]`**

Gives the Levenshtein distance between two strings or two lists: the fewest single-element insertions, deletions and substitutions that turn one into the other. Strings are compared byte by byte, so a multi-byte UTF-8 character counts as several elements.

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

Three edits turn one word into the other

```mathematica
In[5]:= EditDistance["Sunday", "Saturday"]
Out[5]= 3
```

Works on lists, not only strings

```mathematica
In[6]:= EditDistance[{1, 2, 3, 4}, {1, 3, 4}]
Out[6]= 1
```

Three insertions from the empty string

```mathematica
In[7]:= EditDistance["", "abc"]
Out[7]= 3
```

## Implementation notes

**Algorithm.** `builtin_edit_distance` returns the Levenshtein distance — the
fewest single-element insertions, deletions and substitutions turning one
sequence into the other. `dist_levenshtein` runs the standard dynamic program
over two rolling rows rather than the full matrix, so it uses O(min(m, n)) memory
and O(m·n) time. Elements are compared with `expr_eq`, so the one routine serves
both strings (compared character by character) and lists of arbitrary
expressions — `EditDistance[{1, 2, 3}, {1, 3}]` is `1`, as in Mathematica.

**Shape and encoding.** `dist_seq_pair` requires both arguments to be strings or
both to be lists; a mixed pair declines. `dist_seq` either borrows a list's
elements or explodes a string into one `Expr` per *byte* — so a multi-byte UTF-8
character counts as several elements, which matches the ASCII and DNA use cases
and is stated rather than silently assumed.

**Complexity / limits.** O(m·n) time, O(min(m, n)) extra memory; returns `-1`
internally on allocation failure, surfaced as an unevaluated call.
`ATTR_PROTECTED`. See `HammingDistance` for the equal-length positional
analogue.

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

**See also:** [HammingDistance](../../lists-and-iteration/HammingDistance/)

- Source: [`src/list/distance.c`](https://github.com/stblake/mathilda/blob/main/src/list/distance.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)

## Notes & additional examples

### Notes

`EditDistance[a, b]` is the Levenshtein distance: the fewest single-element
insertions, deletions, and substitutions that turn `a` into `b`. Both arguments
must be strings, or both lists — a mixed pair is left unevaluated. Elements are
compared for equality, so the same function serves strings (character by
character) and lists of arbitrary expressions.

Strings are compared byte by byte, so a multi-byte UTF-8 character counts as
several elements. For the positional "how many slots differ" count on
equal-length sequences, use `HammingDistance` instead.
