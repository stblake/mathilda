# StringCount

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringCount["string", "sub"]
    Gives the number of times "sub" appears as a substring of "string".
StringCount["string", patt]
    Gives the number of substrings of "string" matching the string
    expression patt.
StringCount["string", {p1, p2, ...}]
    Counts the occurrences of any of the pi.
StringCount[{s1, s2, ...}, patt]
    Gives the list of results for each of the si.

    Equivalent to Length[StringCases[...]] but does not build the
    matched substrings. Options: Overlaps -> False (default; overlapping
    substrings are not counted as separate), True (overlaps counted
    separately, one substring per start), or All (every matching
    substring at every start); IgnoreCase -> True treats upper/lowercase
    as equivalent.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringCount["the cat sat on the mat", "at"]
Out[1]= 3

In[2]:= StringCount["AAAA", "AA"]
Out[2]= 2

In[3]:= StringCount["AAAA", "AA", Overlaps -> True]
Out[3]= 3

In[4]:= StringCount["AAAA", x__, Overlaps -> All]
Out[4]= 10

In[5]:= StringCount[{"a1", "b22", "ccc"}, DigitCharacter]
Out[5]= {1, 2, 0}

In[6]:= StringCount["aAbB", "a", IgnoreCase -> True]
Out[6]= 2

In[7]:= StringCount["x=1.5, y=-2", NumberString]
Out[7]= 2
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
