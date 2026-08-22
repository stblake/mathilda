# StringReplace

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringReplace["string", patt -> rep]`**

Replaces each non-overlapping match of patt in "string" by rep, with $n replaced by the n-th captured group and $0 by the whole match.

**`StringReplace["string", {patt1 -> rep1, patt2 -> rep2, ...}]`**

Applies a list of replacement rules; at each position the leftmost match wins, ties broken by rule order.

**`StringReplace[{s1, s2, ...}, rules]`**

Gives the list of results for each of the si.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= StringReplace["a13b12c1da32efg", RegularExpression["(\\d+)"] -> "[$1]"]
Out[1]= "a[13]b[12]c[1]da[32]efg"

In[2]:= StringReplace["123 45 6 789", RegularExpression["\\b"] :> "X"]
Out[2]= "X123X X45X X6X X789X"
```

## Algorithm

stringreplace.c - StringReplace[subject, rule | {rules...}]

Replaces every non-overlapping match of a rule's pattern by its right-hand

```text
side, scanning left to right; unmatched text is copied verbatim.  The RHS is
```

a string in which $0/$1... expand to the whole match and capture groups. With several rules, at each position the leftmost match wins, ties broken by

```text
rule order.  A list of subjects threads.
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Characters of 200k chars | 4.38 s | 2.54 s | 0.419 s |
| StringSplit on space, 200k chars | 3.12 s | 4.13 s | 0.723 s |
| StringCases regex, 200k chars | 2.36 s | 3.69 s | 3.34 s |
| StringReplace regex, 200k chars | 1.97 s | 4.74 s | 3.35 s |
| StringReplace literal, 200k chars | 0.332 s | 1.17 s | 0.195 s |
| StringCount substring, 200k chars | 0.224 s | 0.364 s | 0.102 s |

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)
