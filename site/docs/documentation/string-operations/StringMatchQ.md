# StringMatchQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringMatchQ["string", patt]`**

Gives True if the whole "string" matches patt, and False otherwise.

**`StringMatchQ[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si. patt may be RegularExpression\["re"\], a literal string, or a list of alternatives.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= StringMatchQ["12345", RegularExpression["\\d+"]]
Out[1]= True

In[2]:= StringMatchQ[{"12", "x"}, RegularExpression["\\d+"]]
Out[2]= {True, False}
```

### Applications (3)

~~ builds a string expression; ___ is zero or more characters

```mathematica
In[3]:= StringMatchQ["hello", "h" ~~ ___]
Out[3]= True
```

Letters followed by digits

```mathematica
In[4]:= StringMatchQ["abc123", LetterCharacter .. ~~ DigitCharacter ..]
Out[4]= True
```

Threads over the subjects

```mathematica
In[5]:= StringMatchQ[{"cat", "dog"}, "c" ~~ __]
Out[5]= {True, False}
```

## Algorithm

stringmatchq.c - StringMatchQ[subject, pattern]

Returns True if the WHOLE subject string matches the pattern, else False. The pattern may be RegularExpression["re"], a literal string, or a List of

```text
alternatives (matches if any one matches).  A list of subjects threads,
```

giving a list of True/False.

## Implementation notes

**Algorithm.** `builtin_stringmatchq` compiles the pattern *anchored* — `regex_rules_build` with `anchored = 1` wraps each rule as `\A(?:...)\z` — so a match means the pattern covers the whole subject. For each subject it runs `regex_match` at offset 0 and returns `True` on the first rule that matches, else `False`. A `List` pattern becomes several rules (a match if any one matches); a list of subjects threads, with non-string elements copied through unchanged.

The pattern may be a literal string, `RegularExpression[...]`, a general string expression (`~~`, `|`, `..`, character classes, `Pattern`), or a list of alternatives — all handled by the shared translator `wl_pattern_to_regex`.

**Data structures.** A `RegexRule` array, each a compiled `RegexProgram`. Only the whole-match pair is requested (`ov[2]`): a yes/no answer needs no capture pool and no `regex_scan` enumeration.

**Complexity / limits.** One PCRE2 match per rule per subject, early-exiting on the first hit. Byte semantics throughout (no UTF-8 decoding). A non-string subject leaves the call unevaluated.

**Attributes:** `Protected`.

## References

- Source: [`src/strings/regex/stringmatchq.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/stringmatchq.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)

## Notes & additional examples

### Notes

`StringMatchQ` tests whether the *whole* string matches: the pattern is compiled
with a `\A(?:...)\z` wrap so both ends are anchored. The pattern may be a literal
string, `RegularExpression[...]`, a general string expression, or a list of
alternatives (a match if any one matches).

A list of subjects threads, giving a list of `True`/`False`; a non-string subject
leaves the call unevaluated.
