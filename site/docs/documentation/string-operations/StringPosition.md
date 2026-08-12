# StringPosition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringPosition["string", patt]`**

Gives a list of the {start, end} character positions at which substrings matching the string pattern patt occur in "string".

**`StringPosition["string", patt, n]`**

Includes only the first n occurrences.

**`StringPosition["string", {p1, p2, ...}]`**

Gives positions of all the pi.

**`StringPosition[{s1, s2, ...}, patt]`**

Threads over a list of strings. Positions use the form consumed by StringTake / StringReplacePart. Options: Overlaps -\> True (default; overlaps allowed, one substring per start), False (no overlaps), or All (every matching substring); IgnoreCase -\> True treats upper/lowercase as equivalent.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= StringPosition["abXYZaaabXYZaaaaXYZXYZ", "XYZ"]
Out[1]= {{3, 5}, {10, 12}, {17, 19}, {20, 22}}

In[2]:= StringPosition["AABBBAABABBCCCBAAA", x_ ~~ x_]
Out[2]= {{1, 2}, {3, 4}, {4, 5}, {6, 7}, {10, 11}, {12, 13}, {13, 14}, {16, 17}, {17, 18}}
```

### Options (2)

```mathematica
In[3]:= StringPosition["AAAAA", "AA", Overlaps -> False]
Out[3]= {{1, 2}, {3, 4}}

In[4]:= StringPosition["abAB", "a", IgnoreCase -> True]
Out[4]= {{1, 1}, {3, 3}}
```

## Algorithm

stringposition.c - StringPosition[subject, pattern, n]

Returns a List of {start, end} character-position pairs at which substrings of

```text
`subject` match the string pattern `pattern`, in the 1-based inclusive form
consumed by StringTake / StringDrop / StringReplacePart.  The pattern may be a
```

literal string, a general string expression (Blank/Pattern/~~/RegularExpression /character classes), or a List of patterns; a List of subjects threads.

Options:

```text
  Overlaps -> True (default) | False | All
    True  - include overlapping substrings, but only the first (natural) match
            starting at each position.
    False - exclude overlapping substrings (greedy left-to-right, global).
    All   - include every matching substring at every start (all lengths).
  IgnoreCase -> True | False (default)
    Treat upper/lowercase as equivalent.
```

A third positional integer argument n keeps only the first n matches.

The match enumeration itself is regex_scan() in regex_common.c, shared with StringCases and StringCount; this file only turns spans into position pairs.

Byte semantics: like the rest of src/strings, positions are byte offsets (no UTF-8 codepoint decoding), consistent with StringLength / StringPart.

## Implementation notes

**Attributes:** `Protected`.

## See also

[StringTake](../../string-operations/StringTake/), [StringDrop](../../string-operations/StringDrop/), [StringReplacePart](../../string-operations/StringReplacePart/), [StringCases](../../string-operations/StringCases/), [SetOptions](../../assignment-and-rules/SetOptions/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringcount.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcount.c)
- Tests: [`tests/test_stringposition.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringposition.c)
