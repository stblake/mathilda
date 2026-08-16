# StringReplacePart

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringReplacePart["string", "snew", {m, n}]`**

Replaces the characters at positions m through n in "string" by "snew".

**`StringReplacePart["string", "snew", {{m1, n1}, {m2, n2}, ...}]`**

Inserts copies of "snew" at several positions.

**`StringReplacePart["string", {"snew1", "snew2", ...}, {{m1, n1}, ...}]`**

Replaces the characters at each range by the corresponding new string; the two lists must be the same length.

**`StringReplacePart[{s1, s2, ...}, snew, part]`**

Gives the list of results for each of the si.

**`StringReplacePart[new, part]`**

is the operator form: StringReplacePart\[new, part\]\[old\] == StringReplacePart\[old, new, part\]. Positions use the form returned by StringPosition and refer to "string" before any replacement is done. Negative positions count from the end. Positions may not overlap. An empty new string deletes the selected characters.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= StringReplacePart["abcdefghijk", "ABCDEFGH", {2, 5}]
Out[1]= "aABCDEFGHfghijk"

In[2]:= StringReplacePart["abcdefghijk", "ABCDEFGH", {{1, 1}, {3, 5}, {-3, -1}}]
Out[2]= "ABCDEFGHbABCDEFGHfghABCDEFGH"

In[3]:= StringReplacePart["abcdefghijk", "ABCDEFGH", {-3, -2}]
Out[3]= "abcdefghABCDEFGHk"

In[4]:= StringReplacePart["abcdefghijk", {"XYZ", "ABCD"}, {{2, 3}, {-2, -2}}]
Out[4]= "aXYZdefghiABCDk"

In[5]:= StringReplacePart["ABCDEFGH", {2, 5}]["abcdefghijk"]
Out[5]= "aABCDEFGHfghijk"

In[6]:= StringReplacePart["abcde", "", {2, 4}]
Out[6]= "ae"

In[7]:= StringReplacePart["abcde", "XYZ", {{1, 3}, {3, 5}}] StringReplacePart::ovlp: Position {3,5} overlaps previous positions; new string XYZ will not be inserted.

In[8]:= StringReplacePart[] StringReplacePart::argt: StringReplacePart called with 0 arguments; 2 or 3 arguments are expected.
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [StringPosition](../../string-operations/StringPosition/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_strings.c`](https://github.com/stblake/mathilda/blob/main/tests/test_strings.c)
