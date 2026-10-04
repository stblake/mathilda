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

## Examples (7)

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

### Applications (3)

{start, end} pairs, 1-based inclusive

```mathematica
In[5]:= StringPosition["abcabc", "bc"]
Out[5]= {{2, 3}, {5, 6}}
```

Overlaps -> True is the default here

```mathematica
In[6]:= StringPosition["aaaa", "aa"]
Out[6]= {{1, 2}, {2, 3}, {3, 4}}
```

A character-class pattern

```mathematica
In[7]:= StringPosition["a1b2c3", DigitCharacter]
Out[7]= {{2, 2}, {4, 4}, {6, 6}}
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

**Algorithm.** `builtin_stringposition` shares the option-seeding and rule-building machinery of `StringCases`/`StringCount`. `sp_scalar` runs the shared `regex_scan` and turns each span `[ms, me)` into the 1-based inclusive pair `{ms + 1, me}` — the position form that `StringTake`, `StringDrop`, and `StringReplacePart` consume. An optional third integer argument caps the number of pairs returned.

It differs from its two siblings in one default: `Overlaps -> True` (matching the Wolfram Language), so overlapping matches at distinct start positions are listed, where `StringCases`/`StringCount` default to `False`.

**Data structures.** The shared `RegexScan` span array; the result is a `List` of two-element position `List`s. A list of subjects threads.

**Complexity / limits.** Shares the `regex_scan` enumerator, so counts and policies agree with `StringCases`/`StringCount`. Positions are byte offsets (no UTF-8 decoding), consistent with `StringLength`/`StringPart`.

**Attributes:** `Protected`.

## References

**See also:** [StringTake](../../string-operations/StringTake/), [StringDrop](../../string-operations/StringDrop/), [StringReplacePart](../../string-operations/StringReplacePart/), [StringCases](../../string-operations/StringCases/), [SetOptions](../../assignment-and-rules/SetOptions/)

- Source: [`src/strings/regex/stringposition.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/stringposition.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringcount.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcount.c)
- Tests: [`tests/test_stringposition.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringposition.c)

## Notes & additional examples

### Notes

`StringPosition` returns `{start, end}` character positions in the 1-based
inclusive form that `StringTake`, `StringDrop`, and `StringReplacePart` consume.
It shares the `regex_scan` enumerator with `StringCases`/`StringCount`.

Unlike those two it defaults to `Overlaps -> True` (matching the Wolfram
Language), so overlapping matches at distinct starts are listed. An optional
third integer argument keeps only the first `n` matches.
