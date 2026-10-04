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

## Examples (11)

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

### Applications (3)

Replace a character range

```mathematica
In[9]:= StringReplacePart["abcdef", "XY", {2, 3}]
Out[9]= "aXYdef"
```

An empty string deletes the range

```mathematica
In[10]:= StringReplacePart["abcdef", "", {2, 4}]
Out[10]= "aef"
```

One new string per range

```mathematica
In[11]:= StringReplacePart["abcdefgh", {"1", "2"}, {{1, 2}, {5, 6}}]
Out[11]= "1cd2gh"
```

## Implementation notes

**Algorithm.** Two arguments is the operator form `StringReplacePart[new, part][old]`, realised as `Function[StringReplacePart[#1, new, part]]`; three is the direct form (any other arity emits `StringReplacePart::argt`). Each position is a `{m, n}` pair (`srp_parse_range`, negatives → `len + k + 1`), validated against the original string. A single new string is broadcast to every range; a list of new strings must match the range count. Ranges are accepted in order against a `covered` byte-mask — a later range touching an already-claimed byte triggers `StringReplacePart::ovlp` and is dropped — then insertion-sorted by start and assembled in two passes (compute length, then emit literal spans and replacements). A `List` first argument recurses per element.

**Data structures.** An `SrpRange` array (original + resolved positions + borrowed replacement pointer), a `bool[len]` covered-mask, a `const char*` replacement array, and one output buffer.

**Complexity / limits.** `O(nranges · rangelen + len)`. Positions use the `StringPosition` form and refer to the original string; an empty replacement string deletes the selected characters. A malformed/out-of-range position or a new-string/position length mismatch leaves the call unevaluated. Byte-indexed.

**Attributes:** `Protected`.

## References

**See also:** [StringPosition](../../string-operations/StringPosition/)

- Source: [`src/strings/stringreplacepart.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringreplacepart.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_strings.c`](https://github.com/stblake/mathilda/blob/main/tests/test_strings.c)

## Notes & additional examples

### Notes

Position specifications are `{m, n}` first/last character pairs — the form
`StringPosition` returns — with negative positions counting from the end. All
positions refer to the *original* string, before any replacement.

A single new string is broadcast to every range; a list of new strings must match
the number of ranges. Overlapping ranges are not allowed: a later range touching
an already-claimed position triggers `StringReplacePart::ovlp` and is dropped. The
operator form is `StringReplacePart[new, part][old]`.
