# StringPosition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringPosition["string", patt]
    Gives a list of the {start, end} character positions at which
    substrings matching the string pattern patt occur in "string".
StringPosition["string", patt, n]
    Includes only the first n occurrences.
StringPosition["string", {p1, p2, ...}]
    Gives positions of all the pi.
StringPosition[{s1, s2, ...}, patt]
    Threads over a list of strings.

    Positions use the form consumed by StringTake / StringReplacePart.
    Options: Overlaps -> True (default; overlaps allowed, one substring
    per start), False (no overlaps), or All (every matching substring);
    IgnoreCase -> True treats upper/lowercase as equivalent.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringPosition["abXYZaaabXYZaaaaXYZXYZ", "XYZ"]
Out[1]= {{3, 5}, {10, 12}, {17, 19}, {20, 22}}

In[2]:= StringPosition["AABBBAABABBCCCBAAA", x_ ~~ x_]
Out[2]= {{1, 2}, {3, 4}, {4, 5}, {6, 7}, {10, 11}, {12, 13}, {13, 14}, {16, 17}, {17, 18}}

In[3]:= StringPosition["AAAAA", "AA", Overlaps -> False]
Out[3]= {{1, 2}, {3, 4}}

In[4]:= StringPosition["abAB", "a", IgnoreCase -> True]
Out[4]= {{1, 1}, {3, 3}}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
