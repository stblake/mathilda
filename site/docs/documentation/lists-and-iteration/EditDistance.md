# EditDistance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EditDistance[u, v]`**

Gives the Levenshtein distance between two strings or two lists: the fewest single-element insertions, deletions and substitutions that turn one into the other. Strings are compared byte by byte, so a multi-byte UTF-8 character counts as several elements.

## Examples (4)

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

## Implementation notes

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

## See also

[HammingDistance](../../lists-and-iteration/HammingDistance/)

## References

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
