# StringCount

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringCount["string", "sub"]`**

Gives the number of times "sub" appears as a substring of "string".

**`StringCount["string", patt]`**

Gives the number of substrings of "string" matching the string expression patt.

**`StringCount["string", {p1, p2, ...}]`**

Counts the occurrences of any of the pi.

**`StringCount[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si. Equivalent to Length\[StringCases\[...\]\] but does not build the matched substrings. Options: Overlaps -\> False (default; overlapping substrings are not counted as separate), True (overlaps counted separately, one substring per start), or All (every matching substring at every start); IgnoreCase -\> True treats upper/lowercase as equivalent.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringCount["the cat sat on the mat", "at"]
Out[1]= 3

In[2]:= StringCount["AAAA", "AA"]
Out[2]= 2

In[3]:= StringCount[{"a1", "b22", "ccc"}, DigitCharacter]
Out[3]= {1, 2, 0}

In[4]:= StringCount["x=1.5, y=-2", NumberString]
Out[4]= 2
```

### Options (3)

```mathematica
In[5]:= StringCount["AAAA", "AA", Overlaps -> True]
Out[5]= 3

In[6]:= StringCount["AAAA", x__, Overlaps -> All]
Out[6]= 10

In[7]:= StringCount["aAbB", "a", IgnoreCase -> True]
Out[7]= 2
```

## Algorithm

stringcount.c - StringCount[subject, pattern]

```text
Counts the substrings of `subject` that match `pattern`.  This is the
```

counting-only companion of StringCases: it runs the same match enumeration (regex_scan in regex_common.c) but never materialises a substring, so it costs one small span record per match instead of a malloc + Expr + List element.

The pattern may be a literal string, a general string expression (Blank/Pattern/~~/RegularExpression/character classes), a Rule or RuleDelayed (only the LHS matters -- the replacement is irrelevant to a count), or a List

```text
of any of those.  A List of subjects threads, giving one count per subject.
```

Options:

```text
  Overlaps -> False (default) | True | All
    False - overlapping substrings are not counted separately.
    True  - overlapping substrings count separately, but only the first
            matching substring at a given position is counted.
    All   - every matching substring at every position is counted separately.
  IgnoreCase -> True | False (default)
    Treat upper/lowercase as equivalent.
```

Because the scan is shared, StringCount[s, p, opts] is always exactly Length[StringCases[s, p, opts]] and Length[StringPosition[s, p, opts]].

Byte semantics: like the rest of src/strings, offsets are byte offsets (no UTF-8 codepoint decoding), consistent with StringLength / StringPart.

## Implementation notes

**Attributes:** `Protected`.

## See also

[StringCases](../../string-operations/StringCases/), [Rule](../../assignment-and-rules/Rule/), [RuleDelayed](../../assignment-and-rules/RuleDelayed/), [SetOptions](../../assignment-and-rules/SetOptions/), [StringPosition](../../string-operations/StringPosition/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringcount.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcount.c)
