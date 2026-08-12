# RegularExpression

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RegularExpression["regex"]`**

Represents a class of strings given by the PCRE regular expression "regex", for use in StringMatchQ, StringCases, StringReplace and StringSplit.  It is an inert head: it evaluates to itself. Supported syntax includes . \[c1c2\] \[c1-c2\] \[^...\] p\* p+ p? p{m,n}, non-greedy \*? +? ??, groups (...) and alternation |; the classes \d \D \s \S \w \W and \[\[:name:\]\]; the anchors ^ $ \b \B; and inline options (?i) (?m) (?s).  In a replacement right-hand side $n stands for the n-th captured group and $0 for the whole match.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= StringCases["adefgh12c34", RegularExpression["[a-e]+"]]
Out[1]= {"ade", "c"}

In[2]:= StringCases["a23b4222c63333d80", RegularExpression["\\d+"]]
Out[2]= {"23", "4222", "63333", "80"}
```

## Algorithm

regularexpression.c - the RegularExpression[...] head.

RegularExpression["re"] is an inert data head: it carries a PCRE pattern string for use by StringMatchQ / StringCases / StringReplace / StringSplit and does not evaluate to anything else (exactly like Wolfram Language). The builtin therefore validates its argument and otherwise returns NULL so

```text
the expression survives evaluation unchanged.  When the pattern does not
```

compile a RegularExpression::regex diagnostic is emitted (still inert), and when Mathilda was built without PCRE2 a RegularExpression::regavail note is printed the first time an invalid-but-present call is seen.

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

## See also

[StringMatchQ](../../string-operations/StringMatchQ/), [StringCases](../../string-operations/StringCases/), [StringReplace](../../string-operations/StringReplace/), [StringSplit](../../string-operations/StringSplit/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_names.c`](https://github.com/stblake/mathilda/blob/main/tests/test_names.c)
- Tests: [`tests/test_regex.c`](https://github.com/stblake/mathilda/blob/main/tests/test_regex.c)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringcount.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcount.c)
