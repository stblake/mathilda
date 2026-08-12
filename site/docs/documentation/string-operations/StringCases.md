# StringCases

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringCases["string", patt]`**

Gives the list of non-overlapping substrings of "string" that match patt, from left to right.

**`StringCases["string", patt -> rhs]`**

Gives the rhs for each match, with $n replaced by the n-th captured group and $0 by the whole match.

**`StringCases["string", {p1, p2, ...}]`**

Gives the matches of any of the pi.

**`StringCases[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si. Options: Overlaps -\> False (default; overlapping substrings are not treated as separate), True (overlaps separate, one substring per start), or All (every matching substring at every start); IgnoreCase -\> True treats upper/lowercase as equivalent.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= StringCases["a13b12c17a32", RegularExpression["[^a1]"]]
Out[1]= {"3", "b", "2", "c", "7", "3", "2"}

In[2]:= StringCases["AaBBccDDeefG", RegularExpression["[[:upper:]]+"]]
Out[2]= {"A", "BB", "DD", "G"}
```

### Options (3)

```mathematica
In[3]:= StringCases["AAAA", "AA", Overlaps -> True]
Out[3]= {"AA", "AA", "AA"}

In[4]:= StringCases["AAA", x__, Overlaps -> All]
Out[4]= {"AAA", "AA", "A", "AA", "A", "A"}

In[5]:= StringCases["aAbB", "a", IgnoreCase -> True]
Out[5]= {"a", "A"}
```

## Algorithm

stringcases.c - StringCases[subject, pattern]

Returns a List of the substrings of `subject` that match `pattern`, left to

```text
right.  With a rule pattern (patt -> rhs / patt :> rhs) each match is replaced
by the rhs, with $0/$1... expanded to the whole match and capture groups.  The
```

pattern may also be a List of alternatives/rules; at each position the

```text
leftmost match wins, ties broken by rule order.  A list of subjects threads.
```

Options:

```text
  Overlaps -> False (default) | True | All
    False - non-overlapping, greedy left-to-right.
    True  - overlapping substrings count separately, but only the first match
            starting at each position.
    All   - every matching substring at every start (all lengths).
  IgnoreCase -> True | False (default)
    Treat upper/lowercase as equivalent.
```

The match enumeration itself is regex_scan() in regex_common.c, shared with StringCount and StringPosition, so StringCount[s, p, opts] always equals Length[StringCases[s, p, opts]].

## Implementation notes

**Attributes:** `Protected`.

## See also

[SetOptions](../../assignment-and-rules/SetOptions/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringcount.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcount.c)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)
