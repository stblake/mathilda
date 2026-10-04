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

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= StringReplace["a13b12c1da32efg", RegularExpression["(\\d+)"] -> "[$1]"]
Out[1]= "a[13]b[12]c[1]da[32]efg"

In[2]:= StringReplace["123 45 6 789", RegularExpression["\\b"] :> "X"]
Out[2]= "X123X X45X X6X X789X"
```

### Applications (3)

Every match is replaced

```mathematica
In[3]:= StringReplace["hello world", "o" -> "0"]
Out[3]= "hell0 w0rld"
```

A string-pattern rule

```mathematica
In[4]:= StringReplace["a1b2c3", DigitCharacter -> "#"]
Out[4]= "a#b#c#"
```

$n reorders captured groups

```mathematica
In[5]:= StringReplace["2024-01", RegularExpression["(\\d+)-(\\d+)"] -> "$2.$1"]
Out[5]= "01.2024"
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

**Algorithm.** `builtin_stringreplace` builds the rule set *unanchored* and requires every element to carry a replacement RHS (otherwise the call is left unevaluated). `sr_scalar_str` then scans left to right: at each position it takes the leftmost match across all rules (ties broken by rule order), copies the literal text before it, appends the rule's expanded replacement (`regex_rule_replacement`, with `$0`/`$n` substituted), and resumes at the match end. A zero-width match copies one character and advances by one so the loop always makes progress.

**Data structures.** A growable byte buffer (`RegexBuf`); a per-match ovector of up to `REGEX_MAX_PAIRS` pairs holds the whole match and its capture groups for template expansion.

**Complexity / limits.** One left-to-right pass with a PCRE2 probe per rule per position. A list of subjects threads, with non-string elements copied through. Offsets are byte offsets. Only string replacements with `$n` templates are expanded; a non-string RHS falls back to copying the matched text.

**Attributes:** `Protected`.

## References

- Source: [`src/strings/regex/stringreplace.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/stringreplace.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)

## Notes & additional examples

### Notes

`StringReplace` scans left to right, replacing each non-overlapping match by its
rule's right-hand side (a string, with `$0`/`$n` expanded) and copying the
unmatched text. With several rules, the leftmost match wins at each position,
ties broken by rule order.

A zero-width match (e.g. `\b`) inserts the replacement without consuming a
character. A list of subjects threads.
