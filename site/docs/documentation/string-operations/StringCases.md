# StringCases

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringCases["string", patt]
    Gives the list of non-overlapping substrings of "string" that match
    patt, from left to right.
StringCases["string", patt -> rhs]
    Gives the rhs for each match, with $n replaced by the n-th captured
    group and $0 by the whole match.
StringCases[{s1, s2, ...}, patt]
    Gives the list of results for each of the si.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringCases["a13b12c17a32", RegularExpression["[^a1]"]]
Out[1]= {"3", "b", "2", "c", "7", "3", "2"}

In[2]:= StringCases["AaBBccDDeefG", RegularExpression["[[:upper:]]+"]]
Out[2]= {"A", "BB", "DD", "G"}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
