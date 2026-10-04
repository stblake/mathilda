# RegularExpression

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RegularExpression["regex"]`**

Represents a class of strings given by the PCRE regular expression "regex", for use in StringMatchQ, StringCases, StringReplace and StringSplit.  It is an inert head: it evaluates to itself. Supported syntax includes . \[c1c2\] \[c1-c2\] \[^...\] p\* p+ p? p{m,n}, non-greedy \*? +? ??, groups (...) and alternation |; the classes \d \D \s \S \w \W and \[\[:name:\]\]; the anchors ^ $ \b \B; and inline options (?i) (?m) (?s).  In a replacement right-hand side $n stands for the n-th captured group and $0 for the whole match.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= StringCases["adefgh12c34", RegularExpression["[a-e]+"]]
Out[1]= {"ade", "c"}

In[2]:= StringCases["a23b4222c63333d80", RegularExpression["\\d+"]]
Out[2]= {"23", "4222", "63333", "80"}
```

### Applications (4)

A run of digits

```mathematica
In[3]:= StringCases["x1y22z333", RegularExpression["\\d+"]]
Out[3]= {"1", "22", "333"}
```

$n back-references a captured group

```mathematica
In[4]:= StringReplace["John Smith", RegularExpression["(\\w+) (\\w+)"] -> "$2 $1"]
Out[4]= "Smith John"
```

Letters then digits, whole string

```mathematica
In[5]:= StringMatchQ["abc123", RegularExpression["[a-z]+\\d+"]]
Out[5]= True
```

The digits are the delimiters

```mathematica
In[6]:= StringSplit["a1b2c3", RegularExpression["\\d"]]
Out[6]= {"a", "b", "c"}
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

**Algorithm.** `RegularExpression["re"]` is an inert data head: `builtin_regularexpression` validates a single string argument and otherwise returns `NULL`, so the head survives evaluation unchanged and carries raw PCRE source for the consuming builtins (`StringMatchQ`, `StringCases`, `StringReplace`, `StringSplit`, `StringPosition`, …). As a best-effort syntax check it runs `regex_compile` once when PCRE2 is present, emitting `RegularExpression::regex` on a bad pattern while staying inert.

The real work is in the two shared layers the consumers call. The translator `wl_pattern_to_regex` (`src/strings/regex/string_pattern.c`) passes `RegularExpression["re"]` through verbatim as PCRE source (and renders the symbolic string-pattern heads — `Whitespace`, `LetterCharacter`, `~~`, `|`, `..`, `Except`, `Pattern` → capture/backreference — into PCRE too). The engine wrapper `regex_compile` (`src/strings/regex/regex_engine.c`) calls `pcre2_compile` plus a best-effort `pcre2_jit_compile`.

**Data structures.** A `RegexProgram` wraps a `pcre2_code`, a `pcre2_match_data`, and the capture count; the 8-bit PCRE2 library is used because Mathilda strings are byte-oriented. In a replacement RHS, `$0` is the whole match and `$n` the n-th group, expanded by `regex_expand_template`.

**Complexity / limits.** Matching cost is PCRE2's; up to `$0..$63` capture pairs are exposed (`REGEX_MAX_PAIRS`). Backslashes are doubled in a Mathilda string literal (`\\d` → `\d`). Without PCRE2 (`USE_REGEX` unset) the engine stubs out and every string-pattern builtin warns (`regavail`) and stays unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [StringMatchQ](../../string-operations/StringMatchQ/), [StringCases](../../string-operations/StringCases/), [StringReplace](../../string-operations/StringReplace/), [StringSplit](../../string-operations/StringSplit/)

- P. Hazel, *PCRE2 — Perl-Compatible Regular Expressions*, 10.x, https://www.pcre.org/.
- Source: [`src/strings/regex/regularexpression.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/regularexpression.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_names.c`](https://github.com/stblake/mathilda/blob/main/tests/test_names.c)
- Tests: [`tests/test_regex.c`](https://github.com/stblake/mathilda/blob/main/tests/test_regex.c)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringcount.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcount.c)

## Notes & additional examples

### Notes

`RegularExpression["re"]` is an inert head carrying raw PCRE2 source; the string
functions compile it once with `pcre2_compile` and reuse it. In a replacement
right-hand side `$0` is the whole match, `$n` the n-th capture group, and `$$` a
literal `$`.

Backslashes are doubled inside a Mathilda string literal (`\\d` is the single
regex token `\d`). An invalid pattern emits `RegularExpression::regex` but stays
inert; a build without PCRE2 warns and leaves the call unevaluated.
