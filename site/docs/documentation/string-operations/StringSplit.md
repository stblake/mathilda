# StringSplit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringSplit["string"]`**

Splits "string" at runs of whitespace.

**`StringSplit["string", patt]`**

Splits at delimiters matching the string pattern patt.

**`StringSplit["string", {p1, p2, ...}]`**

Splits at any of the pi.

**`StringSplit["string", patt -> val]`**

Inserts val at the position of each delimiter.

**`StringSplit["string", patt, n]`**

Splits into at most n substrings.

**`StringSplit[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si. Empty substrings between adjacent interior delimiters are kept; those at the start or end are dropped unless All is given as the third argument. "" splits at every character. Option: IgnoreCase -\> True.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringSplit["a bbb  cccc aa   d"]
Out[1]= {"a", "bbb", "cccc", "aa", "d"}

In[2]:= StringSplit["a-b:c-d:e-f-g", {":", "-"}]
Out[2]= {"a", "b", "c", "d", "e", "f", "g"}

In[3]:= StringSplit["a b::c d::e f g", "::" -> "--"]
Out[3]= {"a b", "--", "c d", "--", "e f g"}

In[4]:= StringSplit["This is a sentence, which goes on.", Except[WordCharacter] ..]
Out[4]= {"This", "is", "a", "sentence", "which", "goes", "on"}
```

### Applications (3)

Default: split at runs of whitespace

```mathematica
In[5]:= StringSplit["a b c"]
Out[5]= {"a", "b", "c"}
```

A literal delimiter

```mathematica
In[6]:= StringSplit["a,b,c", ","]
Out[6]= {"a", "b", "c"}
```

A character-class delimiter

```mathematica
In[7]:= StringSplit["a1b2c", DigitCharacter]
Out[7]= {"a", "b", "c"}
```

## Algorithm

stringsplit.c - StringSplit[...], the full Wolfram surface.

```text
  StringSplit[s]                 split at runs of whitespace
  StringSplit[s, patt]           split at delimiters matching patt
  StringSplit[s, {p1, p2, ...}]  split at any of the pi
  StringSplit[s, patt -> val]    insert val at each delimiter
  StringSplit[s, patt, n]        at most n substrings
  StringSplit[s, patt, All]      keep leading/trailing empty substrings
  StringSplit[{s1, ...}, patt]   thread over a list of subjects
  IgnoreCase -> True             case-insensitive delimiters
```

The delimiter pattern is translated to PCRE by the shared string-pattern engine (regex_common.c / string_pattern.c), so it accepts literal strings, RegularExpression["re"], the character-class heads (Whitespace, ...), StringExpression (~~), Alternatives (|), Repeated (..), Except, and so on.

Zero-length substrings between two adjacent interior delimiters are kept; empty substrings at the very beginning or end are dropped unless All is given. The empty-string delimiter "" splits at every character.

## Implementation notes

**Algorithm.** `builtin_stringsplit` strips trailing `IgnoreCase` options, accepts 1–3 positional arguments, and compiles the delimiter pattern (default `Whitespace`) to PCRE through the shared engine (`regex_rules_build_ex`). `ss_scalar` scans left to right, pushing each inter-delimiter substring (interior empties kept) and, for a `patt -> val` rule, the inserted value — `make_insertion` expands a string RHS as a `$`-template and binds a named `x:patt :> f[x]` by a scoped `ReplaceAll`. A third argument `n` caps the piece count; `All` keeps the leading/trailing empties that are otherwise trimmed. The empty delimiter `""` short-circuits to a per-byte split (`ss_null`). A list of subjects threads.

**Data structures.** A `PieceVec` of `(Expr, is_sub)` records — the `is_sub` flag marks an empty-trimmable substring versus an inserted rule value — plus the per-match ovector.

**Complexity / limits.** One left-to-right pass with a PCRE2 probe per rule per position; a zero-width delimiter advances by one so the scan progresses. Offsets are byte offsets. The delimiter accepts the full shared string-pattern vocabulary (literals, `RegularExpression`, character classes, `~~`, `|`, `..`, `Except`, …).

**Attributes:** `Protected`.

## References

**See also:** [StringReplace](../../string-operations/StringReplace/)

- Source: [`src/strings/stringsplit.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringsplit.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)

## Notes & additional examples

### Notes

`StringSplit` returns the substrings between non-overlapping matches of the
delimiter, which is translated to PCRE by the shared string-pattern engine
(literals, `RegularExpression`, character classes, `~~`, `|`, `..`, `Except`, …).

Zero-length substrings between two adjacent interior delimiters are kept; leading
and trailing empties are dropped unless `All` is given as a third argument. The
empty-string delimiter `""` splits at every character; `IgnoreCase -> True` folds
case.
